"""Regions BED files, read as Hail intervals for partitioning the combiner and densifying."""

import gzip

import hail as hl
from cpg_utils import to_path


def read_bed_file_as_intervals(bed_path: str) -> list[hl.Interval]:
    """Reads a BED file as GRCh38 intervals.

    Args:
        bed_path: Path to a plain-text or gzipped BED file. Blank lines and '#' lines are skipped.

    Returns:
        One closed interval per BED record, converted from BED's 0-based, half-open coordinates.
    """
    path = to_path(bed_path)
    intervals: list[hl.Interval] = []
    with gzip.open(path.open('rb'), 'rt') if str(path).endswith('.gz') else path.open() as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped or stripped.startswith('#'):
                continue
            chrom, start, end = stripped.split()[:3]
            start_locus = hl.Locus(chrom, int(start) + 1, reference_genome='GRCh38')
            end_locus = hl.Locus(chrom, int(end), reference_genome='GRCh38')
            intervals.append(hl.Interval(start_locus, end_locus, includes_start=True, includes_end=True))
    return intervals
