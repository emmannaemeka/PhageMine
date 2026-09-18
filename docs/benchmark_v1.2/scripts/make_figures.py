"""Regenerate functional yield from validated aggregate evidence (not accuracy)."""
from pathlib import Path
import csv

EXPECTED = {'PhageMine': (804, 439), 'Pharokka': (804, 398), 'Prokka': (675, 328)}


def load_totals(root):
    with (root / 'tables/functional_yield_totals.tsv').open() as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    assert len(rows) == 3 and {r['tool'] for r in rows} == set(EXPECTED)
    for row in rows:
        total, named = EXPECTED[row['tool']]
        assert (int(row['predicted_cds']), int(row['named_product_calls'])) == (total, named)
        assert int(row['hypothetical_or_uncharacterized']) == total - named
        assert abs(float(row['named_fraction']) - named / total) < 0.0000005
    return rows


def main():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root = Path(__file__).resolve().parents[1]
    rows = load_totals(root)
    plt.rcParams.update({'svg.hashsalt': 'phagemine-v12', 'pdf.fonttype': 42})
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    bars = ax.bar([r['tool'] for r in rows], [int(r['named_product_calls']) for r in rows], color=['#0072B2', '#E69F00', '#009E73'])
    ax.bar_label(bars, padding=4)
    ax.set(ylabel='Named product calls', ylim=(0, 510), title='Functional yield — not functional accuracy')
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    for ext in ('png', 'svg', 'pdf'):
        metadata = {'Date': None} if ext == 'svg' else {'CreationDate': None, 'ModDate': None} if ext == 'pdf' else {}
        fig.savefig(root / 'figures' / f'functional_yield.{ext}', dpi=300, metadata=metadata)
    svg = root / 'figures/functional_yield.svg'
    svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
    plt.close(fig)


if __name__ == '__main__':
    main()
