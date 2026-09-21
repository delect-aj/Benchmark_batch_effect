"""
MicroCoRE Framework -- Stage 3
Taxon representation learning from the consensus correlation network

Core idea:
    Treat the consensus correlation matrix as the adjacency matrix of a
    weighted graph (nodes = taxa, edge weights = correlation strengths).
    A graph autoencoder (GAE) learns a low-dimensional embedding for every
    taxon by reconstructing the correlation matrix from inner products of
    embedding vectors.
    An optional phylogenetic smoothness regulariser encourages taxonomically
    close taxa to have similar embeddings.
"""

import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim


class TaxonGraphAutoencoder(nn.Module):
    """
    Minimal graph autoencoder:
        Embedding matrix  E in R^{n_taxa x embed_dim}
        Reconstructed correlations = tanh(E @ E.T)

    tanh is used because correlation coefficients lie in [-1, 1].
    """

    def __init__(self, n_taxa, embed_dim=64):
        super().__init__()
        self.embedding = nn.Parameter(torch.randn(n_taxa, embed_dim) * 0.1)

    def forward(self):
        return torch.tanh(self.embedding @ self.embedding.T)

    def get_embedding(self):
        return self.embedding.detach().cpu().numpy()


def train_graph_autoencoder(consensus_matrix,
                             embed_dim=64,
                             n_epochs=2000,
                             lr=0.01,
                             phylo_matrix=None,
                             smoothness_weight=0.1,
                             device="cpu",
                             print_every=200,
                             random_seed=42):
    """
    Train the graph autoencoder to obtain taxon embeddings.

    Parameters
    ----------
    consensus_matrix : np.ndarray, shape (n_taxa, n_taxa)
        Output of Stage 2.
    embed_dim : int
        Dimensionality of each taxon embedding vector.
    n_epochs : int
    lr : float
        Adam learning rate.
    phylo_matrix : np.ndarray or None
        Optional, shape (n_taxa, n_taxa). Entry (i, j) is the phylogenetic
        proximity between taxa i and j (e.g. 1 = same genus, 0 = unrelated).
        When provided, a Laplacian smoothness term encourages neighbouring
        taxa to have similar embeddings.
    smoothness_weight : float
        Weight of the phylogenetic smoothness regulariser.
    device : str
        'cpu' or 'cuda'.
    print_every : int
    random_seed : int

    Returns
    -------
    taxon_embeddings : np.ndarray, shape (n_taxa, embed_dim)
    final_loss : float
    """
    torch.manual_seed(random_seed)
    np.random.seed(random_seed)

    n_taxa = consensus_matrix.shape[0]
    target = torch.from_numpy(consensus_matrix).float().to(device)

    model = TaxonGraphAutoencoder(n_taxa, embed_dim).to(device)
    optimizer = optim.Adam(model.parameters(), lr=lr)

    # Compute reconstruction loss only on non-zero edges to prevent the
    # large number of zeros from dominating the objective
    edge_mask = (torch.abs(target) > 1e-6).float()
    edge_weight_sum = edge_mask.sum().clamp(min=1.0)

    phylo_tensor = (torch.from_numpy(phylo_matrix).float().to(device)
                    if phylo_matrix is not None else None)

    print(f"  Training graph autoencoder -- "
          f"n_taxa={n_taxa}, embed_dim={embed_dim}, n_epochs={n_epochs}")

    loss_history = []
    for epoch in range(n_epochs):
        model.train()
        optimizer.zero_grad()

        reconstruction = model()

        # Masked MSE reconstruction loss
        recon_loss = ((reconstruction - target) ** 2 * edge_mask).sum() \
                     / edge_weight_sum

        total_loss = recon_loss

        # Phylogenetic smoothness regulariser (Laplacian penalty):
        #   sum_{i,j} A_{ij} * ||e_i - e_j||^2
        if phylo_tensor is not None and smoothness_weight > 0:
            emb = model.embedding
            dist_sq = torch.cdist(emb, emb, p=2) ** 2
            smooth_term = (phylo_tensor * dist_sq).sum() \
                          / phylo_tensor.sum().clamp(min=1.0)
            total_loss = total_loss + smoothness_weight * smooth_term

        total_loss.backward()
        optimizer.step()

        if (epoch + 1) % print_every == 0 or epoch == 0:
            print(f"    Epoch {epoch + 1:5d}/{n_epochs} | "
                  f"recon loss {recon_loss.item():.6f} | "
                  f"total loss {total_loss.item():.6f}")

        loss_history.append(total_loss.item())

    taxon_embeddings = model.get_embedding()
    return taxon_embeddings, loss_history[-1]
