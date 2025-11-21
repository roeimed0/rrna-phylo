# rRNA-Phylo Development Roadmap

## Current Status

### [OK] Completed Features

1. **Phylogenetic Tree Building - FULLY WORKING**
   - UPGMA (distance-based)
   - BioNJ (variance-weighted neighbor-joining)
   - Maximum Likelihood with GTR+Gamma (DNA/RNA)
   - Maximum Likelihood with WAG+Gamma (Protein)
   - See [visualize_all_trees.py](backend/visualize_all_trees.py) for working demo

2. **Distance Calculation**
   - Jukes-Cantor (DNA/RNA)
   - Kimura 2-parameter (DNA/RNA)
   - Poisson correction (Protein)

3. **Sequence Processing**
   - FASTA parsing
   - MUSCLE alignment integration
   - Sequence type detection (DNA/RNA/Protein)
   - Quality validation

4. **ML Optimizations (Level 1-3)**
   - Level 1: Basic GTR model
   - Level 2: GTR + rate heterogeneity
   - Level 3: GTR+Gamma with site pattern compression (CURRENT)

5. **Testing Infrastructure**
   - 11/11 tests passing
   - Comprehensive test coverage
   - [run_all_tests.py](backend/run_all_tests.py) test runner

---

## Remaining Features

### 1. Consensus Tree Methods (HIGH PRIORITY)

**Goal**: Combine multiple trees (UPGMA, BioNJ, ML) into consensus trees with confidence scores.

#### 1.1 Tree Comparison Metrics
- Robinson-Foulds distance (topological difference)
- Branch score distance (weighted by branch lengths)
- Quartet distance
- Path-length distance

**Implementation**: `rrna_phylo/consensus/tree_distance.py`

#### 1.2 Consensus Algorithms
- **Strict Consensus**: Only clades in ALL trees
- **Majority-Rule Consensus**: Clades in >50% of trees
- **Extended Majority-Rule**: Adds compatible clades
- **Greedy Consensus**: Iteratively adds most frequent clades

**Implementation**: `rrna_phylo/consensus/consensus.py`

#### 1.3 Confidence Scoring
- Bootstrap support values
- Frequency-based support (% of trees with each clade)
- Posterior probabilities
- Consistency indices

**Implementation**: `rrna_phylo/consensus/confidence.py`

#### 1.4 Integration
- Update `PhylogeneticTreeBuilder` to include consensus
- Add `trees['consensus']` and `trees['consensus_support']` to output
- Visualize support values on consensus tree

---

### 2. ML Tree Level 4 (HIGH PRIORITY)

**Goal**: Advanced ML optimizations for better accuracy and performance.

#### 2.1 Model Selection
- Automatic model testing (ModelFinder-style)
- AIC/BIC/BICc comparison
- Cross-validation
- Model averaging

**Implementation**: `rrna_phylo/models/model_selection.py`

#### 2.2 Advanced Rate Models
- Branch-specific rate variation
- Free-rate models (FreeRate)
- Site-specific rate categories
- Codon models (if coding sequences)
- Partition models (different rates per region)

**Implementation**: `rrna_phylo/models/ml_tree_level4.py`

#### 2.3 Performance Optimizations
- Vectorized likelihood calculations (NumPy)
- Cython/Numba for inner loops
- Better starting tree heuristics
- Parallel bootstrap analysis
- GPU acceleration (optional, PyTorch/CuPy)

#### 2.4 Tree Search Improvements
- NNI (Nearest Neighbor Interchange)
- SPR (Subtree Pruning and Regrafting)
- TBR (Tree Bisection and Reconnection)
- Stochastic tree perturbation

---

### 3. Gen-AI Integration (MEDIUM PRIORITY)

**Goal**: Use machine learning to enhance tree quality and enable generative capabilities.

#### 3.1 ML-Enhanced Quality Scoring
- Train classifier on tree quality metrics
- Predict phylogenetic signal reliability
- Feature engineering from alignment properties:
  - Sequence similarity distribution
  - Gap patterns
  - Conserved regions
  - Rate variation
- Outlier sequence detection

**Implementation**: `rrna_phylo/ml/quality_scorer.py`

#### 3.2 Multi-Tree Ensemble Methods
- Learn weights for different tree methods
- Weighted consensus based on quality
- Ensemble predictions for branch support
- Tree outlier detection

**Implementation**: `rrna_phylo/ml/ensemble.py`

#### 3.3 Generative Tree Synthesis (ADVANCED)
- Graph Neural Networks (GNN) for tree generation
- Transformers for sequence-to-tree models
- Variational autoencoders for tree space
- Conditional generation with constraints
- Tree completion/refinement

**Implementation**: `rrna_phylo/ml/generative/`

**Dependencies**: PyTorch, DGL/PyTorch Geometric

---

### 4. Backend API (MEDIUM PRIORITY)

**Goal**: REST API for tree building as a service.

#### 4.1 FastAPI Service
- Async request handling
- Job queue (Celery + Redis)
- Progress tracking
- Result caching
- Rate limiting

**Implementation**: `backend/api/main.py`

#### 4.2 API Endpoints

**Tree Building**:
- `POST /api/v1/build-tree` - Submit sequences, get job ID
- `GET /api/v1/job/{job_id}` - Check job status
- `GET /api/v1/results/{job_id}` - Retrieve trees
- `DELETE /api/v1/job/{job_id}` - Cancel/delete job

**Consensus**:
- `POST /api/v1/consensus` - Build consensus from multiple trees
- `POST /api/v1/compare-trees` - Calculate tree distances

**File Upload**:
- `POST /api/v1/upload/fasta` - Upload FASTA file
- `POST /api/v1/upload/tree` - Upload Newick/Nexus tree

#### 4.3 Frontend Integration
- React/Vue web interface (optional)
- Tree visualization (D3.js, Phylo.js)
- Interactive alignment viewer
- Result export (Newick, Nexus, PNG, SVG)

---

## Priority Order

### Phase 1: Core Phylogenetics (Recommended First)
1. **Consensus Tree Methods** - Completes the phylogenetic pipeline
   - Essential for forensics-grade reliability
   - Provides confidence metrics
   - Estimated time: 1-2 weeks

### Phase 2: Advanced ML
2. **ML Tree Level 4** - Enhanced accuracy and performance
   - Model selection framework
   - Better tree search
   - Estimated time: 2-3 weeks

### Phase 3: AI Integration
3. **Gen-AI Integration** - Modern ML enhancements
   - Quality scoring (Phase 3.1)
   - Ensemble methods (Phase 3.2)
   - Generative models (Phase 3.3) - Optional/Research
   - Estimated time: 3-4 weeks (excluding generative)

### Phase 4: Service Deployment
4. **Backend API** - Production deployment
   - FastAPI service
   - Web interface
   - Estimated time: 2-3 weeks

---

## Technical Decisions

### Libraries to Add
- **Consensus/Tree Comparison**: DendroPy (optional) or custom implementation
- **ML Quality Scoring**: scikit-learn, XGBoost
- **Generative Models**: PyTorch, DGL/PyTorch Geometric
- **API**: FastAPI, Celery, Redis
- **Visualization**: matplotlib, seaborn (backend), D3.js (frontend)

### Architecture Considerations
- Keep core phylogenetics pure Python (current approach)
- Add Cython/Numba for performance-critical sections (Level 4)
- ML models as optional add-ons (don't break core functionality)
- API as separate service (microservices pattern)

---

## Next Steps

**Immediate**: Start with Consensus Tree Methods
1. Implement Robinson-Foulds distance
2. Implement majority-rule consensus
3. Add confidence scoring
4. Integrate with `build_trees()` function
5. Add comprehensive tests

**Command to run current demo**:
```bash
cd backend
python visualize_all_trees.py
```

This will show all 9 working tree types (UPGMA/BioNJ/ML for DNA/RNA/Protein).
