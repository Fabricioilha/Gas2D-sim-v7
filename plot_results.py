#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gera gráficos científicos a partir de traj_cm.csv:
- |v_CM|(t) em m/s
- T_eff(t) em K

Uso:
    python3 plot_results.py
    python3 plot_results.py --csv traj_cm.csv --out out_plots --smooth-sec 0.5

Requisitos: matplotlib (sem seaborn, sem pandas).
"""

import argparse
import csv
import math
import os
import sys
from typing import List, Tuple

import matplotlib.pyplot as plt


def read_csv_cm(path: str) -> Tuple[List[float], List[float], List[float]]:
    """
    Lê CSV no formato: t,vel_cm,Teff (com cabeçalho).
    Retorna listas: t [s], vcm [m/s], Teff [K].
    """
    t, vcm, teff = [], [], []
    with open(path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header is None:
            raise ValueError("CSV vazio.")
        # Aceita tanto nomes esperados quanto posição das colunas
        # Esperado: t, vel_cm, Teff
        # Fallback: tenta converter as 3 primeiras colunas como floats
        idx_t = 0
        idx_v = 1
        idx_T = 2 if len(header) >= 3 else None

        # Detecta por nome (case-insensitive, ignora espaços)
        norm = [h.strip().lower() for h in header]
        try:
            idx_t = norm.index("t")
        except ValueError:
            pass
        for cand in ("vel_cm", "|v_cm|", "vcm"):
            if cand in norm:
                idx_v = norm.index(cand)
                break
        for cand in ("teff", "t_eff", "temperature", "temperature_eff", "temp_eff"):
            if cand in norm:
                idx_T = norm.index(cand)
                break

        for row in reader:
            if not row or all((c is None or c.strip() == "") for c in row):
                continue
            try:
                tt = float(row[idx_t])
                vv = float(row[idx_v])
                if idx_T is not None and idx_T < len(row):
                    TT = float(row[idx_T])
                else:
                    TT = float("nan")
            except (ValueError, IndexError):
                # pula linhas malformadas
                continue
            t.append(tt)
            vcm.append(vv)
            teff.append(TT)
    if len(t) < 2:
        raise ValueError("Poucos pontos no CSV para plotar.")
    return t, vcm, teff


def moving_average(y: List[float], win_samples: int) -> List[float]:
    """
    Média móvel simples com janela centrada (tamanho ímpar preferível).
    Para bordas, usa janela truncada (sem padding).
    """
    if win_samples <= 1:
        return y[:]
    n = len(y)
    out = [0.0] * n
    half = win_samples // 2
    csum = [0.0] * (n + 1)
    for i in range(n):
        csum[i + 1] = csum[i] + (0.0 if math.isnan(y[i]) else y[i])
    for i in range(n):
        a = max(0, i - half)
        b = min(n - 1, i + half)
        cnt = (b - a + 1)
        s = csum[b + 1] - csum[a]
        out[i] = s / max(1, cnt)
    return out


def infer_dt(t: List[float]) -> float:
    """Estima dt médio ignorando outliers simples."""
    if len(t) < 2:
        return float("nan")
    diffs = [t[i + 1] - t[i] for i in range(len(t) - 1)]
    diffs = [d for d in diffs if d > 0]
    if not diffs:
        return float("nan")
    diffs.sort()
    # mediana robusta
    mid = len(diffs) // 2
    if len(diffs) % 2 == 1:
        return diffs[mid]
    return 0.5 * (diffs[mid - 1] + diffs[mid])


def ensure_outdir(outdir: str) -> None:
    if outdir and not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)


def scientific_style(ax, xlabel: str, ylabel: str, title: str):
    """Configuração de estilo científico limpa (sem cores fixas)."""
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, which="both", linestyle="--", linewidth=0.6, alpha=0.6)
    ax.tick_params(axis="both", which="major", labelsize=10)
    # limites automáticos; sem legendas poluídas por padrão


def save_both(fig: plt.Figure, outdir: str, basename: str):
    png_path = os.path.join(outdir, f"{basename}.png") if outdir else f"{basename}.png"
    pdf_path = os.path.join(outdir, f"{basename}.pdf") if outdir else f"{basename}.pdf"
    fig.tight_layout()
    fig.savefig(png_path, dpi=200, bbox_inches="tight")
    fig.savefig(pdf_path, bbox_inches="tight")
    print(f"Salvo: {png_path}")
    print(f"Salvo: {pdf_path}")


def plot_vcm(t: List[float], vcm: List[float], outdir: str, smooth_win: int = 1):
    yy = vcm
    if smooth_win > 1:
        yy = moving_average(vcm, smooth_win)

    fig = plt.figure(figsize=(7, 4.2))
    ax = fig.add_subplot(111)
    ax.plot(t, yy, linewidth=1.5)
    scientific_style(ax,
                     xlabel="Tempo, t (s)",
                     ylabel=r"$|\mathbf{v}_{CM}(t)|$ (m s$^{-1}$)",
                     title=r"Módulo da velocidade do centro de massa")
    save_both(fig, outdir, "vcm_vs_t")
    plt.close(fig)


def plot_teff(t: List[float], teff: List[float], outdir: str, smooth_win: int = 1):
    # Pode haver NaN se o CSV antigo não tinha Teff
    yy = teff
    if all(math.isnan(v) for v in teff):
        print("Aviso: coluna Teff ausente ou inválida no CSV; gráfico de T_eff será pulado.")
        return
    if smooth_win > 1:
        # substitui NaN por média local simples antes de suavizar
        base = [0.0 if math.isnan(v) else v for v in teff]
        yy = moving_average(base, smooth_win)

    fig = plt.figure(figsize=(7, 4.2))
    ax = fig.add_subplot(111)
    ax.plot(t, yy, linewidth=1.5)
    scientific_style(ax,
                     xlabel="Tempo, t (s)",
                     ylabel=r"$T_{\mathrm{eff}}(t)$ (K)",
                     title=r"Temperatura efetiva (2D)")
    save_both(fig, outdir, "teff_vs_t")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description="Plot de |v_CM|(t) e T_eff(t) a partir de traj_cm.csv")
    parser.add_argument("--csv", type=str, default="traj_cm.csv",
                        help="caminho do CSV (padrão: traj_cm.csv)")
    parser.add_argument("--out", type=str, default="plots",
                        help="pasta de saída (padrão: plots)")
    parser.add_argument("--smooth-sec", type=float, default=0.0,
                        help="suavização por média móvel (janela em segundos). 0 desativa.")
    parser.add_argument("--max-points", type=int, default=0,
                        help="opcional: máximo de pontos (downsample uniforme). 0 = sem limite.")
    args = parser.parse_args()

    if not os.path.isfile(args.csv):
        print(f"Erro: não encontrei '{args.csv}'. Execute primeiro o simulador para gerar o CSV.")
        sys.exit(1)

    try:
        t, vcm, teff = read_csv_cm(args.csv)
    except Exception as e:
        print(f"Erro ao ler CSV: {e}")
        sys.exit(1)

    # Downsample opcional para arquivos enormes (mantém forma global)
    if args.max_points and args.max_points > 0 and len(t) > args.max_points:
        step = max(1, len(t) // args.max_points)
        t = t[::step]
        vcm = vcm[::step]
        teff = teff[::step]

    # Suavização: converte janela de segundos para amostras
    dt = infer_dt(t)
    if args.smooth_sec > 0 and math.isfinite(dt) and dt > 0:
        win_samples = max(1, int(round(args.smooth_sec / dt)))
        if win_samples % 2 == 0:  # prefira janela ímpar
            win_samples += 1
    else:
        win_samples = 1

    ensure_outdir(args.out)
    plot_vcm(t, vcm, args.out, smooth_win=win_samples)
    plot_teff(t, teff, args.out, smooth_win=win_samples)

    print("Concluído com sucesso.")


if __name__ == "__main__":
    main()
