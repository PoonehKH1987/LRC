# LRC: Logistic Regression Model and Clustering Approach

## Overview
LRC (Logistic Regression Model and Clustering Approach) is a MATLAB-based computational framework designed for predicting conformational B-cell epitopes. The method integrates physicochemical, structural, and statistical properties with logistic regression and clustering to enhance epitope prediction accuracy. 

### Workflow:
1. **Antigen Surface Graph Construction**: The antigen surface is modeled as a weighted graph using `Antigen_Surface_Graph.m`, where vertices represent surface residues and edges define their spatial relationships.
2. **Feature Extraction & Weight Assignment**: Residues are characterized based on physicochemical (hydrophobicity, polarity, flexibility, antigenicity), structural (protrusion index), and statistical properties. Logistic regression (`LR.m`) is used to compute residue weights.
3. **Graph-Based Clustering**: The Markov Cluster Algorithm (`MCL.m`, `MCLCluster.m`) is applied to segment the weighted graph. Clusters with lower average weights are filtered out (`Filtering.m`).
4. **Epitope Prediction**: Residues in the highest-weighted clusters are identified as potential epitope residues using `LRC.m`.

This approach is validated using a dataset of 90 Antibody-Antigen complexes, with training and testing sets carefully selected to ensure robust performance evaluation.

## Features
- Integrates antigen surface modeling, feature extraction, and clustering for epitope prediction
- Implements logistic regression for residue weighting
- Uses the Markov Cluster Algorithm (MCL) for graph-based clustering
- Designed for bioinformatics applications, particularly immunoinformatics and vaccine design

## Prerequisites
Before using this repository, ensure you have the following:
- MATLAB R2020a or later
- Statistics and Machine Learning Toolbox (for logistic regression and clustering functions)
- Optimization Toolbox (optional, for model tuning)

## Installation
1. Clone the repository to your local machine:
   ```bash
   git clone https://github.com/PoonehKH1987/LRC.git
   ```
2. Open MATLAB and navigate to the project directory:
   ```matlab
   cd 'path_to_LRC'
   ```
3. Ensure all necessary toolboxes are installed by running:
   ```matlab
   ver
   ```

## Usage
### 1. Construct Antigen Surface Graph
Run:
```matlab
Protein = Antigen_Surface_Graph4('3NH7.pdb', 'A');
```

### 2. Compute Feature Weights
```matlab
[W, W_graph] = weigthed_Gragh(Protein.Matrix, X);
```

### 3. Apply Markov Clustering
```matlab
Clusters = MCL(W_graph, Protein.Vertices, 30);
```

### 4. Predict Epitope Residues
```matlab
E = LRC('GetArea-ASA3NH7.txt', '3NH7.pdb', 'A', X);
```

### 5. Evaluate Model Performance
```matlab
accuracy = sum(E == Y_test) / length(Y_test);
fprintf('Model Accuracy: %.2f%%\n', accuracy * 100);
```

## Example Applications
LRC can be applied to:
- **Immunoinformatics**: Predicting B-cell epitopes for vaccine development
- **Medical Diagnosis**: Identifying antigenic regions for disease diagnostics
- **Bioinformatics**: Analyzing protein-protein interactions and immune responses

## Contact
For questions or suggestions, contact Pooneh Khodabakhsh at poonehkhodabakhsh87@gmail.com.
