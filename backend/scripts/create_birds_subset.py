#!/usr/bin/env python3
"""
Create a curated birds-only subset of Arcosauria test dataset
Keeps 40 diverse bird species + 1 reptile outgroup (turtle)
"""

from Bio import SeqIO

# Read all sequences
input_file = "c:/Users/User/Desktop/projects/rrna-phylo/backend/data/test/Arcosauria_test.fasta"
output_file = "c:/Users/User/Desktop/projects/rrna-phylo/backend/data/test/birds_test.fasta"

all_seqs = list(SeqIO.parse(input_file, "fasta"))

# Reptiles (non-bird) - we'll keep just ONE as outgroup
reptiles = [
    "Acontias_percivali",      # Lizard
    "Mauremys_sinensis",       # Chinese pond turtle
    "Trachemys_scripta"        # Red-eared slider turtle
]

# Select diverse bird families (40 species representing major groups)
selected_birds = [
    # Ratites (flightless)
    "Struthio_camelus",        # Ostrich
    "Rhea_americana",          # Rhea
    "Dromaius_novaehollandiae", # Emu
    "Apteryx_australis",       # Kiwi

    # Galliformes (gamebirds)
    "Gallus_gallus",           # Chicken
    "Meleagris_gallopavo",     # Turkey
    "Coturnix_coturnix",       # Quail

    # Anseriformes (waterfowl)
    "Anas_platyrhynchos",      # Mallard duck

    # Columbiformes (pigeons/doves)
    "Columba_livia",           # Rock pigeon

    # Accipitriformes (raptors)
    "Haliaeetus_leucocephalus", # Bald eagle
    "Falco_cherrug",           # Saker falcon
    "Neophron_percnopterus",   # Egyptian vulture

    # Strigiformes (owls)
    "Bubo_virginianus",        # Great horned owl

    # Piciformes (woodpeckers)
    "Picoides_pubescens",      # Downy woodpecker

    # Psittaciformes (parrots)
    "Melopsittacus_undulatus", # Budgerigar

    # Apodiformes (swifts/hummingbirds)
    "Apus_apus",               # Common swift
    "Anthracothorax_nigricollis", # Black-throated mango

    # Coraciiformes (kingfishers/hornbills)
    "Coracias_caudatus",       # Lilac-breasted roller
    "Upupa_epops",             # Hoopoe

    # Passeriformes (perching birds - most diverse!)
    "Passer_montanus",         # Eurasian tree sparrow
    "Corvus_corax",            # Common raven (if available) or Pica pica
    "Pica_pica",               # Magpie
    "Turdus_merula",           # Blackbird (if available) or similar thrush
    "Catharus_guttatus",       # Hermit thrush
    "Sturnus_vulgaris",        # Starling (if available) or similar
    "Hirundo_aethiopica",      # Ethiopian swallow
    "Phylloscopus_trochiloides", # Greenish warbler
    "Sylvia_sp",               # Warbler
    "Cinclus_pallasii",        # Brown dipper
    "Lanius_isabellinus",      # Isabelline shrike
    "Bombycilla_cedrorum",     # Cedar waxwing

    # Gruiformes (cranes/rails)
    "Grus_canadensis",         # Sandhill crane

    # Charadriiformes (shorebirds/gulls)
    "Charadrius_semipalmatus", # Semipalmated plover
    "Larus_glaucoides",        # Iceland gull

    # Pelecaniformes (pelicans/herons)
    "Pelecanus_occidentalis",  # Brown pelican

    # Suliformes (boobies/cormorants)
    "Sula_sula",               # Red-footed booby
    "Phalacrocorax_varius",    # Pied cormorant
]

# Find sequences
kept_sequences = []
kept_names = []

# Add turtle outgroup first
for seq in all_seqs:
    if "Trachemys_scripta" in seq.description:
        kept_sequences.append(seq)
        kept_names.append("Trachemys_scripta (OUTGROUP)")
        break

# Add selected birds
for bird in selected_birds:
    for seq in all_seqs:
        if bird in seq.description and bird not in kept_names:
            kept_sequences.append(seq)
            kept_names.append(bird)
            break

print(f"Selected {len(kept_sequences)} sequences:")
print(f"  - 1 reptile outgroup: {kept_names[0]}")
print(f"  - {len(kept_sequences)-1} birds")
print()

# Write output
SeqIO.write(kept_sequences, output_file, "fasta")
print(f"✅ Created: {output_file}")
print(f"   Total sequences: {len(kept_sequences)}")
