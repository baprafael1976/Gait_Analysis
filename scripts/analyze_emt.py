"""
analyze_emt.py
==============
Análise de forças na marcha a partir de arquivo .emt com duas plataformas de força.

Uso:
    python analyze_emt.py <arquivo.emt> [--report outputs/report.csv]
                                        [--fz-threshold 10]
                                        [--save-figures]

Colunas esperadas no arquivo (separadas por espaço):
    Frame | Time | Fx1 | Fy1 | Fz1 | Fx2 | Fy2 | Fz2

Unidades: força em N, tempo em s.
Frequência de amostragem assumida: 500 Hz.
"""

import argparse
import os
import math

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import signal

# ---------------------------------------------------------------------------
# Parâmetros globais
# ---------------------------------------------------------------------------
CUT_FREQ = 10           # Hz — frequência de corte do filtro passa-baixa
SAMPLE_FREQ = 500       # Hz — frequência de amostragem
FORCE_PLATE_LENGTH = 0.6  # m — comprimento de cada plataforma
NUM_FORCE_PLATES = 4    # número de plataformas percorridas

REPORT_KEYS = {
    "total_time":       "Total Time [s]",
    "number_of_steps":  "Number of Steps",
    "sym_idx_fx":       "Symmetry Index Fx",
    "sym_idx_fy":       "Symmetry Index Fy",
    "sym_idx_fz":       "Symmetry Index Fz",
    "integral_r":       "Integral Fy (R) [N·s]",
    "integral_l":       "Integral Fy (L) [N·s]",
    "mech_energ_r":     "Mechanical Energy Expenditure (R) [J]",
    "mech_energ_l":     "Mechanical Energy Expenditure (L) [J]",
    "single_support_r": "Single Support Phase (R) [s]",
    "single_support_l": "Single Support Phase (L) [s]",
    "double_support":   "Double Support Phase [s]",
    "velocity_ms":      "Average Velocity [m/s]",
    "velocity_kmh":     "Average Velocity [km/h]",
    "step_freq":        "Step Frequency [step/s]",
    "step_time_num":    "Step Time (Number of Steps Based) [s]",
    "step_time_sup":    "Step Time (Support Phase Based) [s]",
}

# ---------------------------------------------------------------------------
# Funções auxiliares
# ---------------------------------------------------------------------------

def remove_zero(values):
    """Remove zeros de uma lista numérica."""
    return [v for v in values if v != 0]


def get_average(values):
    """Retorna a média de uma lista; retorna 0.0 se a lista estiver vazia."""
    if not values:
        return 0.0
    return sum(values) / len(values)


# ---------------------------------------------------------------------------
# Filtro Butterworth passa-baixa
# ---------------------------------------------------------------------------

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = signal.butter(order, normal_cutoff, btype='low', analog=False)
    return b, a


def butter_lowpass_filter(data, cutoff=CUT_FREQ, fs=SAMPLE_FREQ, order=5):
    """
    Aplica filtro de Butterworth passa-baixa segmento a segmento,
    respeitando regiões de NaN (ausência de contato).
    Segmentos com menos de 21 amostras são deixados sem filtragem.
    """
    b, a = butter_lowpass(cutoff, fs, order=order)
    arr = np.array(data, dtype=float)
    isnan = np.isnan(arr)

    if not isnan.any():
        return signal.filtfilt(b, a, arr)

    filtered = arr.copy()
    edges = np.where(np.diff(isnan.astype(int)) != 0)[0] + 1
    edges = np.concatenate([[0], edges, [len(arr)]])

    for start, end in zip(edges[:-1], edges[1:]):
        segment = arr[start:end]
        if not np.isnan(segment).any() and (end - start) > 20:
            filtered[start:end] = signal.filtfilt(b, a, segment)

    return filtered


# ---------------------------------------------------------------------------
# Carregamento de dados
# ---------------------------------------------------------------------------

def load_data(file_path, fz_threshold=0.0):
    """
    Lê o arquivo .emt e retorna um DataFrame filtrado.

    Linhas onde Fz1 e Fz2 ficam abaixo de `fz_threshold` (sem contato em
    ambas as plataformas simultaneamente) são descartadas.  Quando
    `fz_threshold` é 0, o comportamento original (baseado em NaN) é mantido.
    """
    rows = []
    with open(file_path, 'r') as fh:
        for line in fh:
            values = line.split()
            if len(values) != 8:
                continue
            rows.append(values)

    df = pd.DataFrame(rows, columns=['Frame', 'Time', 'Fx1', 'Fy1', 'Fz1', 'Fx2', 'Fy2', 'Fz2'])
    df = df.apply(pd.to_numeric, errors='coerce')

    if fz_threshold > 0:
        # Descarta linhas onde ambas as plataformas estão sem contato
        both_off = (df['Fz1'].fillna(0).abs() < fz_threshold) & \
                   (df['Fz2'].fillna(0).abs() < fz_threshold)
        df = df[~both_off].reset_index(drop=True)
    else:
        # Comportamento original: descarta linhas onde Fz1 e Fz2 são NaN
        both_nan = df['Fz1'].isna() & df['Fz2'].isna()
        df = df[~both_nan].reset_index(drop=True)

    for col in ['Fx1', 'Fy1', 'Fz1', 'Fx2', 'Fy2', 'Fz2']:
        df[col] = butter_lowpass_filter(df[col].values)

    return df


# ---------------------------------------------------------------------------
# Gráficos
# ---------------------------------------------------------------------------

def plot_raw_data(raw_data, file_path, save_dir=None):
    """Plota as três componentes de força ao longo do tempo."""
    fig, axes = plt.subplots(3, 1, figsize=(15, 10))
    label = os.path.basename(file_path)

    components = [
        ('Fx1', 'Fx2', 'Força em X', 'Force X (N)'),
        ('Fy1', 'Fy2', 'Força em Y', 'Force Y (N)'),
        ('Fz1', 'Fz2', 'Força em Z', 'Force Z (N)'),
    ]
    for ax, (col1, col2, title, ylabel) in zip(axes, components):
        ax.plot(raw_data['Time'], raw_data[col1], label=f'{col1} (Right Foot)')
        ax.plot(raw_data['Time'], raw_data[col2], label=f'{col2} (Left Foot)')
        ax.set_title(f'{title} — {label}')
        ax.set_xlabel('Time (s)')
        ax.set_ylabel(ylabel)
        ax.legend()

    fig.tight_layout()
    _save_or_show(fig, save_dir, label, 'all_forces')


def plot_y_direction(raw_data, file_path, save_dir=None):
    """Plota apenas a componente Y da força."""
    label = os.path.basename(file_path)
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(raw_data['Time'], raw_data['Fy1'], label='Fy1 (Right Foot)')
    ax.plot(raw_data['Time'], raw_data['Fy2'], label='Fy2 (Left Foot)')
    ax.set_title(f'Force in Y Direction — {label}')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Force Y (N)')
    ax.legend()
    fig.tight_layout()
    _save_or_show(fig, save_dir, label, 'force_y')


def _save_or_show(fig, save_dir, file_label, suffix):
    """Salva a figura em `save_dir` ou exibe na tela."""
    if save_dir:
        os.makedirs(save_dir, exist_ok=True)
        base = os.path.splitext(file_label)[0]
        path = os.path.join(save_dir, f'{base}_{suffix}.png')
        fig.savefig(path, dpi=150)
        print(f"Figura salva em: {path}")
        plt.close(fig)
    else:
        plt.show()


# ---------------------------------------------------------------------------
# Cálculos biomecânicos
# ---------------------------------------------------------------------------

def get_gait_symmetry(raw_data):
    """Índices de simetria por pico para Fx, Fy e Fz."""
    peaks = {col: raw_data[col].max() for col in ['Fx1', 'Fy1', 'Fz1', 'Fx2', 'Fy2', 'Fz2']}

    def sym_idx(a, b):
        denom = peaks[a] + peaks[b]
        return (2 * abs(peaks[a] - peaks[b]) / denom) if denom != 0 else float('nan')

    return {
        'sym_idx_fx': sym_idx('Fx1', 'Fx2'),
        'sym_idx_fy': sym_idx('Fy1', 'Fy2'),
        'sym_idx_fz': sym_idx('Fz1', 'Fz2'),
    }


def get_integrals(raw_data):
    """Integral de Fy para cada pé via regra do trapézio."""
    valid = raw_data.dropna(subset=['Fy1', 'Fy2'])
    integral_r = np.trapz(valid['Fy1'], valid['Time'])
    integral_l = np.trapz(valid['Fy2'], valid['Time'])
    return integral_r, integral_l


def get_mechanical_energy(raw_data):
    """
    Estimativa de gasto mecânico a partir de integrais força-tempo.
    Interprete como estimativa grosseira: não modela massa nem cinemática real.
    """
    valid = raw_data.dropna()
    if valid.empty:
        return 0.0, 0.0

    dt = np.diff(valid['Time'].values)
    dt = np.insert(dt, 0, dt[0])

    total_energy = []
    for fx, fy, fz in [('Fx1', 'Fy1', 'Fz1'), ('Fx2', 'Fy2', 'Fz2')]:
        vx = np.cumsum(valid[fx].values * dt)
        vy = np.cumsum(valid[fy].values * dt)
        vz = np.cumsum(valid[fz].values * dt)
        dx = np.cumsum(vx * dt)
        dy = np.cumsum(vy * dt)
        dz = np.cumsum(vz * dt)
        t = valid['Time'].values
        work = (np.trapz(valid[fx].values * dx, t) +
                np.trapz(valid[fy].values * dy, t) +
                np.trapz(valid[fz].values * dz, t))
        total_energy.append(work)

    return total_energy[0], total_energy[1]


def get_support_phase(raw_data, r_col='Fy1', l_col='Fy2'):
    """
    Estima as fases de suporte simples (direito/esquerdo) e duplo a partir
    de NaN nas colunas de força vertical.

    Retorna três listas de durações em segundos: direito, esquerdo, duplo.
    """
    right_phases, left_phases, both_phases = [], [], []
    r_count = l_count = b_count = 0
    prev_state = None

    for _, row in raw_data.iterrows():
        r_nan = math.isnan(float(row[r_col])) if not pd.isna(row[r_col]) else True
        l_nan = math.isnan(float(row[l_col])) if not pd.isna(row[l_col]) else True

        if not r_nan and l_nan:
            state = 'right'
        elif r_nan and not l_nan:
            state = 'left'
        elif not r_nan and not l_nan:
            state = 'both'
        else:
            state = 'none'

        if state != prev_state:
            if state == 'right':
                right_phases.append(0)
            elif state == 'left':
                left_phases.append(0)
            elif state == 'both':
                both_phases.append(0)
            prev_state = state

        if state == 'right' and right_phases:
            right_phases[-1] += 1
        elif state == 'left' and left_phases:
            left_phases[-1] += 1
        elif state == 'both' and both_phases:
            both_phases[-1] += 1

    to_seconds = lambda lst: [v / SAMPLE_FREQ for v in remove_zero(lst)]
    return to_seconds(right_phases), to_seconds(left_phases), to_seconds(both_phases)


def get_number_of_steps(single_right, single_left):
    return len(single_right) + len(single_left)


def get_step_frequency(number_of_steps, total_time):
    return number_of_steps / total_time if total_time > 0 else 0.0


def get_step_time_num(total_time, number_of_steps):
    return total_time / number_of_steps if number_of_steps > 0 else 0.0


def get_step_time_sup(single_right, single_left, double):
    avg_r = get_average(single_right)
    avg_l = get_average(single_left)
    avg_d = get_average(double)
    return (avg_r + avg_l + avg_d) / 2


def get_velocity(total_time):
    distance = NUM_FORCE_PLATES * FORCE_PLATE_LENGTH
    return distance / total_time if total_time > 0 else 0.0


# ---------------------------------------------------------------------------
# Relatório
# ---------------------------------------------------------------------------

def export_report(data, report_path):
    """Acrescenta uma linha ao CSV de relatório (cria o arquivo se necessário)."""
    os.makedirs(os.path.dirname(report_path) or '.', exist_ok=True)
    header = not os.path.isfile(report_path)
    data.to_csv(report_path, mode='a', header=header, index=False)


def print_report(data):
    print("\n--- Relatório ---")
    for key, label in REPORT_KEYS.items():
        if key in data.columns:
            val = data[key].iloc[0]
            print(f"  {label}: {float(val):.4f}")
    print()


# ---------------------------------------------------------------------------
# Análise principal
# ---------------------------------------------------------------------------

def analyze_file(file_path, report_path='outputs/report.csv',
                 fz_threshold=0.0, save_figures=False):
    """
    Lê, filtra e analisa um arquivo .emt de plataforma de força.
    Gera gráficos, imprime resultados e exporta um CSV de relatório.
    """
    save_dir = 'outputs' if save_figures else None

    raw_data = load_data(file_path, fz_threshold=fz_threshold)

    plot_raw_data(raw_data, file_path, save_dir=save_dir)
    plot_y_direction(raw_data, file_path, save_dir=save_dir)

    time_clean = raw_data.dropna(subset=['Time'])
    total_time = float(time_clean['Time'].iloc[-1] - time_clean['Time'].iloc[0])

    symmetry = get_gait_symmetry(raw_data)
    integral_r, integral_l = get_integrals(raw_data)
    mech_r, mech_l = get_mechanical_energy(raw_data)

    single_r, single_l, double = get_support_phase(raw_data)
    n_steps = get_number_of_steps(single_r, single_l)
    step_freq = get_step_frequency(n_steps, total_time)
    step_time_num = get_step_time_num(total_time, n_steps)
    step_time_sup = get_step_time_sup(single_r, single_l, double)
    velocity_ms = get_velocity(total_time)

    report_data = pd.DataFrame([{
        'name':             os.path.basename(file_path),
        'total_time':       total_time,
        'number_of_steps':  n_steps,
        **symmetry,
        'integral_r':       integral_r,
        'integral_l':       integral_l,
        'mech_energ_r':     mech_r,
        'mech_energ_l':     mech_l,
        'single_support_r': get_average(single_r),
        'single_support_l': get_average(single_l),
        'double_support':   get_average(double),
        'velocity_ms':      velocity_ms,
        'velocity_kmh':     velocity_ms * 3.6,
        'step_freq':        step_freq,
        'step_time_num':    step_time_num,
        'step_time_sup':    step_time_sup,
    }])

    export_report(report_data.drop('name', axis=1), report_path)
    print_report(report_data)


# ---------------------------------------------------------------------------
# Ponto de entrada
# ---------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description='Análise de forças na marcha a partir de arquivo .emt')
    parser.add_argument('file_path', help='Caminho para o arquivo .emt')
    parser.add_argument('--report', default='outputs/report.csv',
                        help='Caminho do CSV de relatório (padrão: outputs/report.csv)')
    parser.add_argument('--fz-threshold', type=float, default=0.0,
                        help='Limiar de Fz (N) para detecção de contato. '
                             '0 = usa NaN do arquivo (padrão)')
    parser.add_argument('--save-figures', action='store_true',
                        help='Salva os gráficos em outputs/ em vez de exibir na tela')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    analyze_file(
        file_path=args.file_path,
        report_path=args.report,
        fz_threshold=args.fz_threshold,
        save_figures=args.save_figures,
    )
