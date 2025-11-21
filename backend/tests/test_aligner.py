import tempfile
import os
from rrna_phylo.io.aligner import MuscleAligner
from rrna_phylo.io.fasta_parser import Sequence

# Determine project root relative to current file
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
LOCAL_MUSCLE = os.path.join(PROJECT_ROOT, "muscle.exe")


def create_unaligned_sequences():
    """Create test sequences in a temporary FASTA file."""
    content = """>seq1 Escherichia_coli_16S
ATGCATGCATGCATGC
>seq2 Bacillus_subtilis_16S
ATGCTGCGATGCATGCTT
>seq3 Staphylococcus_aureus_16S
ATGCATGCGATGC
>seq4 Pseudomonas_aeruginosa_16S
ATGATGCATGCATGCATG
"""
    temp = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta", mode="w")
    temp.write(content)
    temp.close()
    return temp.name


def get_muscle_path():
    """Return path to local muscle.exe if exists, else None to use default."""
    if os.path.exists(LOCAL_MUSCLE):
        return LOCAL_MUSCLE
    return None  # Let MuscleAligner use its default path


def test_muscle_installation():
    print("\n" + "=" * 60)
    print("TEST: MUSCLE Installation & Auto-Install")
    print("=" * 60)

    aligner = MuscleAligner(muscle_executable=get_muscle_path())

    installed = aligner._check_muscle_installed()
    version = getattr(aligner, "version_string", "N/A")

    print(f"\nMUSCLE installed: {installed}")
    print(f"MUSCLE version:   {version}\n")
    return installed


def test_alignment():
    print("\n" + "=" * 60)
    print("TEST: MSA with MUSCLE")
    print("=" * 60)

    input_file = create_unaligned_sequences()
    tmp = tempfile.NamedTemporaryFile(delete=False)
    output_file = tmp.name
    tmp.close()

    try:
        aligner = MuscleAligner(muscle_executable=get_muscle_path())
        aligned_sequences = aligner.align(input_file, output_file)

        print("\n[OK] Alignment completed")
        for s in aligned_sequences:
            print(f"  {s.id}: {s.aligned_length} bp aligned")

        lengths = [s.aligned_length for s in aligned_sequences]
        print(f"\nAligned length uniform: {len(set(lengths)) == 1}")
        print(f"Aligned length:         {lengths[0]}")

        return True

    except Exception as e:
        print(f"\n✗ Alignment failed: {e}")
        return False

    finally:
        if os.path.exists(input_file):
            os.remove(input_file)
        if os.path.exists(output_file):
            os.remove(output_file)


def test_align_sequence_objects():
    print("\n" + "=" * 60)
    print("TEST: Align Sequence Objects")
    print("=" * 60)

    sequences = [
        Sequence("s1", "Test1", "ATGCATGCATGC"),
        Sequence("s2", "Test2", "ATGCTGCGATGCTT"),
        Sequence("s3", "Test3", "ATGCGATGC"),
    ]

    output_file = tempfile.NamedTemporaryFile(delete=False).name

    try:
        aligner = MuscleAligner(muscle_executable=get_muscle_path())
        aligned_sequences = aligner.align_sequences(sequences, output_file)

        print(f"\n[OK] Aligned {len(aligned_sequences)} sequences")
        for seq in aligned_sequences:
            print(f"  {seq.id}: aligned length = {seq.aligned_length}")

        return True

    except Exception as e:
        print(f"\n✗ Failed: {e}")
        return False

    finally:
        if os.path.exists(output_file):
            os.remove(output_file)


def main():
    print("\n" + "=" * 60)
    print("RUNNING MUSCLE TEST SUITE")
    print("=" * 60)

    test_muscle_installation()
    test_alignment()
    test_align_sequence_objects()

    print("\n" + "=" * 60)
    print("ALL TESTS COMPLETE [OK]")
    print("=" * 60)


if __name__ == "__main__":
    main()
