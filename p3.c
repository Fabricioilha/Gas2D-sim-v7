#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

// ---------------------- Constantes REVISADAS ----------------------
#define KB 1.380649e-23                 // J/K
#define T_BASE 333.15                   // Temperatura da base em Kelvin (60 graus Celsius)
#define ALPHA 0.0001                    // Fator adimensional
#define Q_PER_HIT (ALPHA * KB * T_BASE) // Energia adicionada por colisão com a base
#define DT 1e-4                         // Passo de tempo (microsegundos)
#define TOTAL_TIME 100                  // Tempo total de simulação em segundos
#define LARGURA 0.10                    // Largura da caixa em metros (10 cm)
#define ALTURA 0.10                     // Altura da caixa em metros (10 cm)
#define MASS 4.65e-26                   // Massa realista de partículas (aproximadamente)
#define N_PARTICULAS 300                // Número de partículas
#define SEED 42                         // Semente RNG
#define M_PI 3.14159265358979323846     // Definição de PI

// ---------------------- Estruturas -------------------------------
typedef struct {
    double x, y;                        // posição
    double vx, vy;                      // velocidade
} Particle;

typedef struct {
    double x, y;                        // posição do CM
    double vx, vy;                      // velocidade do CM
    double temp_media;                  // temperatura média do sistema
    double energia_total;               // energia total do sistema
} CenterOfMass;

// ---------------------- Utilidades --------------------------------
double urand() { return (double)rand() / (double)RAND_MAX; }
double urand_range(double a, double b) { return a + (b - a) * urand(); }

// Inicializa partículas com velocidades fixas e direções aleatórias
void init_particles(Particle *p, int N) {
    double velocidade_fixa = 0.1;  // 10 cm/s
    
    for (int i = 0; i < N; ++i) {
        p[i].x = urand_range(0.0, LARGURA);
        p[i].y = urand_range(0.0, ALTURA);
        
        // Velocidade fixa mas direção aleatória
        double angle = urand_range(0, 2 * M_PI);
        p[i].vx = velocidade_fixa * cos(angle);
        p[i].vy = velocidade_fixa * sin(angle);
    }
}

// Colisões com paredes e base
void verificar_colisao(Particle *p) {
    const double EPSX = 1e-12;
    const double EPSY = 1e-12;

    // --- Laterais (x): elásticas
    if (p->x <= 0.0) {
        p->x = EPSX;
        p->vx = fabs(p->vx);           // rebate para dentro
    } else if (p->x >= LARGURA) {
        p->x = LARGURA - EPSX;
        p->vx = -fabs(p->vx);          // rebate para dentro
    }

    // --- Base (y=0): elástica + aquecimento
    if (p->y <= 0.0) {
        p->y = EPSY;

        // energia antes/depois
        double v2      = p->vx*p->vx + p->vy*p->vy;
        double E_atual = 0.5 * MASS * v2;
        double E_nova  = E_atual + Q_PER_HIT;
        if (E_nova < 0.0) E_nova = 0.0;

        if (E_atual > 1e-30) {
            // A) aquecimento isotrópico (escala vx e vy)
            double fator = sqrt(E_nova / E_atual);
            p->vx *= fator;
            p->vy *= fator;
            if (p->vy <= 0.0) p->vy = -p->vy;  // garante que saia para cima
        } else {
            // sem energia anterior: sorteia estado saindo para cima (±60°)
            double v_nova = sqrt(2.0 * E_nova / MASS);
            double angle  = urand_range(-M_PI/3.0, M_PI/3.0);
            p->vx = v_nova * sin(angle);
            p->vy = v_nova * cos(angle);       // > 0
        }
    }

    // --- Teto (y = ALTURA): elástica
    if (p->y >= ALTURA) {
        p->y  = ALTURA - EPSY;
        p->vy = -fabs(p->vy);
    }
}

static inline void step_position(Particle *p) {
    p->x += p->vx * DT;
    p->y += p->vy * DT;
}

CenterOfMass calcular_cm(const Particle *p, int N) {
    CenterOfMass cm = {0,0,0,0,0,0};
    double soma_v2 = 0.0;
    double soma_energia = 0.0;
    
    for (int i = 0; i < N; ++i) {
        cm.x  += p[i].x;
        cm.y  += p[i].y;
        cm.vx += p[i].vx;
        cm.vy += p[i].vy;
        
        double v2 = p[i].vx * p[i].vx + p[i].vy * p[i].vy;
        soma_v2 += v2;
        soma_energia += 0.5 * MASS * v2;
    }
    
    cm.x  /= N;
    cm.y  /= N;
    cm.vx /= N;
    cm.vy /= N;
    
    // Calcula temperatura média (para 2 graus de liberdade)
    cm.temp_media = (MASS * soma_v2) / (2.0 * N * KB);
    cm.energia_total = soma_energia;
    
    return cm;
}

int main() {
    srand(SEED);
    
    // Calcula algumas quantidades para debug
    double v_thermal = sqrt(2.0 * KB * T_BASE / MASS);
    
    printf("=== SIMULAÇÃO DE GÁS COM BASE QUENTE ===\n");
    printf("N = %d partículas\n", N_PARTICULAS);
    printf("T_BASE = %.1f K (%.1f°C)\n", T_BASE, T_BASE - 273.15);
    printf("Massa por partícula = %.2e kg\n", MASS);
    printf("Q_PER_HIT = %.2e J\n", Q_PER_HIT);
    printf("Velocidade térmica estimada = %.2e m/s\n", v_thermal);
    printf("Tempo total: %d s | DT: %.1e s\n", TOTAL_TIME, DT);
    printf("Dimensões da caixa: %.2f x %.2f m\n\n", LARGURA, ALTURA);

    Particle *part = (Particle*)malloc(N_PARTICULAS * sizeof(Particle));
    if (!part) { perror("malloc"); return 1; }

    init_particles(part, N_PARTICULAS);
    
    FILE *fout = fopen("dados_simulacao.csv", "w");
    if (!fout) { perror("dados_simulacao.csv"); free(part); return 1; }
    
    fprintf(fout, "tempo, xCM, yCM, vxCM, vyCM, vCM, temperatura, energia\n");

    int steps = (int)lrint(TOTAL_TIME / DT);
    int report_interval = steps / 10;
    
    printf("Progresso da simulação: ");
    fflush(stdout);
    
    for (int step = 0; step <= steps; ++step) {
        double t = step * DT;

        if (report_interval > 0 && step % report_interval == 0) {
            printf("#");
            fflush(stdout);
        }

        // Atualiza posições
        for (int i = 0; i < N_PARTICULAS; ++i) 
            step_position(&part[i]);
        
        // Verifica colisões
        for (int i = 0; i < N_PARTICULAS; ++i) 
            verificar_colisao(&part[i]);

        // Calcula métricas
        CenterOfMass cm = calcular_cm(part, N_PARTICULAS);
        double vcm = sqrt(cm.vx*cm.vx + cm.vy*cm.vy);
        
        fprintf(fout, "%.6f, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e, %.6e\n", 
                t, cm.x, cm.y, cm.vx, cm.vy, vcm, cm.temp_media, cm.energia_total);
    }
    
    printf(" Concluído!\n");
    fclose(fout);
    free(part);

    // Gera scripts gnuplot - CORREÇÃO: Layout 3x1 (vertical)
    FILE *gnu = fopen("plot_tudo.gnu", "w");
    if (gnu) {
        fprintf(gnu, "set datafile separator ','\n");
        fprintf(gnu, "set terminal png size 1600,2400\n");  // Mais alto para 3 gráficos verticais
        fprintf(gnu, "set output 'analise_completa.png'\n");
        fprintf(gnu, "set multiplot layout 3,1\n");  // 3 linhas, 1 coluna
        fprintf(gnu, "set grid\n\n");
        
        // Gráfico 1: Velocidade
        fprintf(gnu, "set title 'Velocidade do Centro de Massa' font ',14'\n");
        fprintf(gnu, "set xlabel 'Tempo (s)' font ',12'\n");
        fprintf(gnu, "set ylabel '|v| (m/s)' font ',12'\n");
        fprintf(gnu, "set yrange [0:*]\n");
        fprintf(gnu, "plot 'dados_simulacao.csv' using 1:6 with lines linewidth 3 title 'Módulo da Velocidade'\n\n");
        
        // Gráfico 2: Temperatura
        fprintf(gnu, "set title 'Temperatura do Sistema' font ',14'\n");
        fprintf(gnu, "set xlabel 'Tempo (s)' font ',12'\n");
        fprintf(gnu, "set ylabel 'Temperatura (K)' font ',12'\n");
        fprintf(gnu, "set yrange [0:*]\n");
        fprintf(gnu, "plot 'dados_simulacao.csv' using 1:7 with lines linewidth 3 title 'Temperatura Média', \\\n");
        fprintf(gnu, "%f with lines linewidth 2 linetype -1 title 'T_{base} = %.1f K'\n\n", T_BASE, T_BASE);
        
        // Gráfico 3: Energia
        fprintf(gnu, "set title 'Energia Total do Sistema' font ',14'\n");
        fprintf(gnu, "set xlabel 'Tempo (s)' font ',12'\n");
        fprintf(gnu, "set ylabel 'Energia (J)' font ',12'\n");
        fprintf(gnu, "set format y '%%g'\n");
        fprintf(gnu, "set yrange [0:*]\n");
        fprintf(gnu, "plot 'dados_simulacao.csv' using 1:8 with lines linewidth 3 title 'Energia Total'\n\n");
        
        fprintf(gnu, "unset multiplot\n");
        fclose(gnu);
    }

    // Script apenas para velocidade - CORREÇÃO: Gráfico largo
    FILE *gnu_vel = fopen("plot_velocidade.gnu", "w");
    if (gnu_vel) {
        fprintf(gnu_vel, "set datafile separator ','\n");
        fprintf(gnu_vel, "set terminal png size 1600,800\n");  // Largo e baixo
        fprintf(gnu_vel, "set output 'velocidade_cm.png'\n");
        fprintf(gnu_vel, "set title 'Evolução da Velocidade do Centro de Massa\\nN=%d partículas, T_{base}=%.1f K, α=%.1e' font ',14'\n", 
                N_PARTICULAS, T_BASE, ALPHA);
        fprintf(gnu_vel, "set xlabel 'Tempo (s)' font ',12'\n");
        fprintf(gnu_vel, "set ylabel '|v| (m/s)' font ',12'\n");
        fprintf(gnu_vel, "set grid\n");
        fprintf(gnu_vel, "set key top right\n");
        fprintf(gnu_vel, "set yrange [0:*]\n");
        fprintf(gnu_vel, "plot 'dados_simulacao.csv' using 1:6 with lines linewidth 3 title 'Módulo da Velocidade do Centro de Massa'\n");
        fclose(gnu_vel);
    }

    // Script para cada gráfico individual - CORREÇÃO: Largura total
    FILE *gnu_individual = fopen("plot_individual.gnu", "w");
    if (gnu_individual) {
        fprintf(gnu_individual, "set datafile separator ','\n");
        fprintf(gnu_individual, "set terminal png size 1600,800\n");
        
        // Gráfico de velocidade individual
        fprintf(gnu_individual, "set output 'velocidade_individual.png'\n");
        fprintf(gnu_individual, "set title 'Velocidade do Centro de Massa\\nN=%d partículas, T_{base}=%.1f K' font ',16'\n", 
                N_PARTICULAS, T_BASE);
        fprintf(gnu_individual, "set xlabel 'Tempo (s)' font ',14'\n");
        fprintf(gnu_individual, "set ylabel '|v| (m/s)' font ',14'\n");
        fprintf(gnu_individual, "set grid linewidth 1\n");
        fprintf(gnu_individual, "set key top right font ',12'\n");
        fprintf(gnu_individual, "set yrange [0:*]\n");
        fprintf(gnu_individual, "plot 'dados_simulacao.csv' using 1:6 with lines linewidth 4 title 'Velocidade do CM'\n\n");
        
        // Gráfico de temperatura individual
        fprintf(gnu_individual, "set output 'temperatura_individual.png'\n");
        fprintf(gnu_individual, "set title 'Temperatura do Sistema\\nN=%d partículas, T_{base}=%.1f K' font ',16'\n", 
                N_PARTICULAS, T_BASE);
        fprintf(gnu_individual, "set xlabel 'Tempo (s)' font ',14'\n");
        fprintf(gnu_individual, "set ylabel 'Temperatura (K)' font ',14'\n");
        fprintf(gnu_individual, "set grid linewidth 1\n");
        fprintf(gnu_individual, "set key top right font ',12'\n");
        fprintf(gnu_individual, "set yrange [0:*]\n");
        fprintf(gnu_individual, "plot 'dados_simulacao.csv' using 1:7 with lines linewidth 4 title 'Temperatura Média', \\\n");
        fprintf(gnu_individual, "%f with lines linewidth 2 linetype -1 title 'Temperatura da Base'\n\n", T_BASE);
        
        // Gráfico de energia individual
        fprintf(gnu_individual, "set output 'energia_individual.png'\n");
        fprintf(gnu_individual, "set title 'Energia Total do Sistema\\nN=%d partículas' font ',16'\n", N_PARTICULAS);
        fprintf(gnu_individual, "set xlabel 'Tempo (s)' font ',14'\n");
        fprintf(gnu_individual, "set ylabel 'Energia (J)' font ',14'\n");
        fprintf(gnu_individual, "set format y '%%g'\n");
        fprintf(gnu_individual, "set grid linewidth 1\n");
        fprintf(gnu_individual, "set key top right font ',12'\n");
        fprintf(gnu_individual, "set yrange [0:*]\n");
        fprintf(gnu_individual, "plot 'dados_simulacao.csv' using 1:8 with lines linewidth 4 title 'Energia Total'\n");
        
        fclose(gnu_individual);
    }

    printf("\n=== ARQUIVOS GERADOS ===\n");
    printf("dados_simulacao.csv   - Dados da simulação\n");
    printf("plot_tudo.gnu         - Script para 3 gráficos verticais\n");
    printf("plot_velocidade.gnu   - Script apenas para velocidade (largo)\n");
    printf("plot_individual.gnu   - Script para gráficos individuais\n");
    printf("\nExecute:\n");
    printf("  gnuplot plot_tudo.gnu          # 3 gráficos verticais\n");
    printf("  gnuplot plot_velocidade.gnu    # Apenas velocidade\n");
    printf("  gnuplot plot_individual.gnu    # Gráficos individuais separados\n");
    
    return 0;
}