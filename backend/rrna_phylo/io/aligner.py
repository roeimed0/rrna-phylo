import os
import subprocess
from rrna_phylo.io.fasta_parser import Sequence, FastaParser
from typing import List, Optional
import tempfile

class MuscleAligner:
    """MUSCLE wrapper compatible with v3 and v5."""

    def __init__(self, muscle_executable: Optional[str] = None):
        # Default to project root muscle.exe
        # __file__ is in rrna_phylo/io/, so go up 3 levels: io -> rrna_phylo -> backend -> project_root
        project_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
        default_path = os.path.join(project_root, "muscle.exe")
        self.muscle_executable = muscle_executable or default_path

        if not os.path.exists(self.muscle_executable):
            raise FileNotFoundError(
                f"MUSCLE not found at '{self.muscle_executable}'. "
                "Place muscle.exe in the project root."
            )

        # Detect version
        r = subprocess.run([self.muscle_executable, "-version"], capture_output=True, text=True)
        self.version_string = r.stdout.strip()
        try:
            self.version_major = int(self.version_string.split()[1].split(".")[0])
        except Exception:
            self.version_major = 3  # fallback to v3 syntax

    def _check_muscle_installed(self) -> bool:
        try:
            r = subprocess.run([self.muscle_executable, "-version"], capture_output=True, text=True)
            return r.returncode == 0
        except Exception:
            return False

    def align(self, input_fasta: str, output_fasta: str, max_iterations: Optional[int] = None) -> List[Sequence]:
        """Align sequences in a FASTA file."""
        # v5 uses -align / -output
        if self.version_major >= 5:
            cmd = [self.muscle_executable, "-align", input_fasta, "-output", output_fasta]
        else:
            cmd = [self.muscle_executable, "-in", input_fasta, "-out", output_fasta]

        if max_iterations is not None:
            cmd.extend(["-maxiters", str(max_iterations)])

        # Timeout: 30 minutes (1800s) for large datasets (100+ sequences)
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=1800)
        if r.returncode != 0:
            raise RuntimeError(f"MUSCLE alignment failed:\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}")

        parser = FastaParser()
        return parser.parse(output_fasta)

    def align_sequences(self, sequences: List[Sequence], output_file: str) -> List[Sequence]:
        """Align a list of Sequence objects in their original order."""
        with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".fasta") as tmp:
            for seq in sequences:
                tmp.write(f">{seq.id} {seq.description}\n{seq.sequence}\n")
            tmp_path = tmp.name

        aligned = self.align(tmp_path, output_file)
        os.remove(tmp_path)

        # CRITICAL FIX: MUSCLE reorders sequences alphabetically!
        # We must restore original input order to prevent taxon mismatches
        # Create mapping from ID to aligned sequence
        id_to_aligned_seq = {seq.id: seq for seq in aligned}

        # Re-sort to match original input order
        aligned_sorted = []
        for original_seq in sequences:
            if original_seq.id in id_to_aligned_seq:
                aligned_sorted.append(id_to_aligned_seq[original_seq.id])
            else:
                raise ValueError(f"Sequence {original_seq.id} not found in aligned output!")

        # Preserve unique display names from original sequences
        # Create a mapping from id to unique_display_name
        id_to_unique_name = {}
        for seq in sequences:
            if hasattr(seq, '_unique_display_name') and seq._unique_display_name:
                id_to_unique_name[seq.id] = seq._unique_display_name

        # Apply unique names to aligned sequences
        for aligned_seq in aligned_sorted:
            if aligned_seq.id in id_to_unique_name:
                aligned_seq._unique_display_name = id_to_unique_name[aligned_seq.id]

        return aligned_sorted
