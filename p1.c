#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
// gcc p1.c -o p1 -lm

/* ================== PARÂMETROS ================== */
#define DT               0.01            /* s  (use pequeno p/ estabilidade do Euler) */
#define TEMPO_SIMULACAO  600.0             /* s  */
#define NUM_PARTICULAS   10000             /* N  */
#define M_PI 3.14159265358979323846

const double LX = 1.0;                    /* m  largura da caixa */
const double LY = 1.0;                    /* m  altura da caixa  */

/* Físicas */
const double KB = 1.380649e-23;           /* J/K */
const double MASSA = 4.65e-26;            /* kg (≈ N2; ajuste à vontade) */
const double T_BASE = 278.15;             /* K  (parede isotérmica na base y=0) */

/* Saída */
const char* ARQ_CM   = "traj_cm.csv";         /* t, |v_CM|, T_eff */
const char* ARQ_SNAP = "traj_particulas.csv"; /* snapshots (opcional, arquivo grande) */

/* ================== ESTRUTURAS ================== */
typedef struct {
    double x, y;   /* posição atual (m) */
    double vx, vy; /* velocidade atual (m/s) */
    double m;      /* massa (kg) */
} Particula;

/* ================== UTILITÁRIOS ================== */
static inline double tempo_colisao_eixo(double pos, double vel, double L) {
    const double EPS = 1e-15;
    if (vel >  EPS) return (L - pos) / vel;   /* indo p/ parede máxima */
    if (vel < -EPS) return (0.0 - pos) / vel; /* indo p/ parede zero   */
    return INFINITY;                           /* sem movimento no eixo */
}

/* Parede térmica fixa (y=0): define v_y = sqrt(k_B T / m) sempre para cima.
   - vx é preservado (parede sem atrito tangencial).
   - Sem aleatoriedade: comportamento determinístico.
*/
static void aplicar_parede_termica_base(Particula* p) {
    const double vy_thermal = sqrt(KB * T_BASE / p->m);  /* m/s */
    p->vy = vy_thermal;  /* garantidamente para cima */
    /* p->vx inalterado */
}

/* Avança UMA partícula por Euler, cortando o passo nos choques e aplicando física de parede */
static void step_euler(Particula* p, double dt) {
    double restante = dt;
    const double EPS = 1e-15;

    while (restante > EPS) {
        /* tempo até próxima colisão (x e y) */
        double tx = tempo_colisao_eixo(p->x, p->vx, LX);
        double ty = tempo_colisao_eixo(p->y, p->vy, LY);
        double tcol = fmin(tx, ty);

        if (!isfinite(tcol) || tcol < 0.0 || tcol > restante) {
            /* Euler puro no subpasso inteiro */
            p->x += p->vx * restante;
            p->y += p->vy * restante;

            /* saneamento numérico: manter no domínio */
            if (p->x < 0) p->x = 0; else if (p->x > LX) p->x = LX;
            if (p->y < 0) p->y = 0; else if (p->y > LY) p->y = LY;

            restante = 0.0;
        } else {
            /* avança até o impacto */
            if (tcol > EPS) {
                p->x += p->vx * tcol;
                p->y += p->vy * tcol;
            }

            /* quais paredes? (empate possível) */
            int bateu_x = fabs(tcol - tx) < 1e-14;
            int bateu_y = fabs(tcol - ty) < 1e-14;

            /* grampeia exatamente na parede */
            if (bateu_x) p->x = (p->vx > 0 ? LX : 0.0);
            if (bateu_y) p->y = (p->vy > 0 ? LY : 0.0);

            /* colisão elástica: inverte componente normal */
            if (bateu_x) p->vx = -p->vx;
            if (bateu_y) {
                p->vy = -p->vy;

                /* se foi na BASE (y==0), aplica a parede térmica fixa */
                if (p->y == 0.0) {
                    aplicar_parede_termica_base(p);
                }
            }

            restante -= tcol;
        }
    }
}

/* Inicializa N partículas "frias"; a base aquecerá o sistema com o tempo */
static void inicializar(Particula* ps, int n, unsigned seed) {
    srand(seed);
    const double V0 = 0.05; /* m/s */

    for (int i = 0; i < n; i++) {
        ps[i].x = ((double)rand()/RAND_MAX) * LX;
        ps[i].y = ((double)rand()/RAND_MAX) * LY;

        double ang = 2.0 * M_PI * ((double)rand()/RAND_MAX);
        ps[i].vx = V0 * cos(ang);
        ps[i].vy = V0 * sin(ang);

        ps[i].m = MASSA;
    }
}

/* Centro de massa instantâneo (pos e vel) */
static void cm(const Particula* ps, int n, double* x, double* y, double* vx, double* vy) {
    double mt=0, sx=0, sy=0, svx=0, svy=0;
    for (int i=0;i<n;i++){
        mt  += ps[i].m;
        sx  += ps[i].m * ps[i].x;
        sy  += ps[i].m * ps[i].y;
        svx += ps[i].m * ps[i].vx;
        svy += ps[i].m * ps[i].vy;
    }
    if (mt>0){
        *x = sx/mt; *y = sy/mt; *vx = svx/mt; *vy = svy/mt;
    } else {
        *x = *y = *vx = *vy = 0.0;
    }
}

/* Temperatura efetiva 2D (equipartição): (1/NkB) * sum_i (1/2 m_i v_i^2) */
static double temperatura_efetiva(const Particula* ps, int n) {
    if (n <= 0) return 0.0;
    double sum_m_v2 = 0.0;
    for (int i=0;i<n;i++){
        double v2 = ps[i].vx*ps[i].vx + ps[i].vy*ps[i].vy;
        sum_m_v2 += ps[i].m * v2;
    }
    return (sum_m_v2) / (2.0 * KB * (double)n);
}

int main(void) {
    Particula* ps = (Particula*) malloc(NUM_PARTICULAS * sizeof(Particula));
    if (!ps) { fprintf(stderr, "Falha ao alocar partículas\n"); return 1; }

    inicializar(ps, NUM_PARTICULAS, 1234u);

    FILE* fcm = fopen(ARQ_CM, "w");
    if (!fcm) { fprintf(stderr, "Erro ao abrir %s\n", ARQ_CM); free(ps); return 1; }
    fprintf(fcm, "t,vel_cm,Teff\n");  /* cabeçalho atualizado */

    FILE* fsnap = fopen(ARQ_SNAP, "w"); /* opcional (pode ser grande) */
    if (fsnap) {
        fprintf(fsnap, "t");
        for (int i=0;i<NUM_PARTICULAS;i++)
            fprintf(fsnap, ",x_%d,y_%d,vx_%d,vy_%d", i,i,i,i);
        fprintf(fsnap, "\n");
    }

    int passos = (int)floor(TEMPO_SIMULACAO/DT) + 1;
    for (int k=0;k<passos;k++){
        double t = k*DT;

        for (int i=0;i<NUM_PARTICULAS;i++)
            step_euler(&ps[i], DT);

        /* |v_CM| */
        double xcm,ycm,vxcm,vycm;
        cm(ps, NUM_PARTICULAS, &xcm,&ycm,&vxcm,&vycm);
        double vcm = hypot(vxcm, vycm);

        /* T_eff */
        double Teff = temperatura_efetiva(ps, NUM_PARTICULAS);

        /* grava saída compacta */
        fprintf(fcm, "%.9f,%.9e,%.9e\n", t, vcm, Teff);

        /* snapshots esparsos */
        if (fsnap && (k % 50 == 0)) {
            fprintf(fsnap, "%.9f", t);
            for (int i=0;i<NUM_PARTICULAS;i++)
                fprintf(fsnap, ",%.9f,%.9f,%.9e,%.9e", ps[i].x, ps[i].y, ps[i].vx, ps[i].vy);
            fprintf(fsnap, "\n");
        }
    }

    if (fsnap) fclose(fsnap);
    fclose(fcm);
    free(ps);

    printf("OK! Salvo '%s' (t, |v_CM|, T_eff). Snapshots em '%s'.\n", ARQ_CM, ARQ_SNAP);
    return 0;
}
