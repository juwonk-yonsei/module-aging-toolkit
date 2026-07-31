"""Attention-MIL core model + training utilities (batched/masked, GPU-efficient).

Bag = (donor) set of instance feature vectors (instance = sub-pseudobulk module
scores). Weakly-supervised binary bag classification.

Two pooling modes with IDENTICAL capacity elsewhere, so any difference is
attributable to the aggregation:
  - "attention": gated attention pooling (Ilse et al., 2018)
  - "mean":      plain (masked) mean pooling

Bags of variable size are zero-padded into a (B, maxN, D) batch with a boolean
mask, so a whole dataset is one batched forward per epoch.

GPU: physical GPU 2 (RTX 4090) pinned via PCI-bus ordering.
"""
import os
os.environ.setdefault("CUDA_DEVICE_ORDER", "PCI_BUS_ID")
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "2")
import numpy as np
import torch
import torch.nn as nn

DEVICE = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
NEG_INF = -1e9


class MIL(nn.Module):
    def __init__(self, in_dim, hidden=64, att_dim=32, pool="attention", dropout=0.25):
        super().__init__()
        self.pool = pool
        self.encoder = nn.Sequential(
            nn.Linear(in_dim, hidden), nn.ReLU(), nn.Dropout(dropout))
        if pool == "attention":
            self.att_V = nn.Linear(hidden, att_dim)
            self.att_U = nn.Linear(hidden, att_dim)
            self.att_w = nn.Linear(att_dim, 1)
        self.head = nn.Linear(hidden, 1)

    def forward(self, X, mask):
        """X: (B, N, in_dim); mask: (B, N) bool (True = valid instance).
        Returns logits (B,) and attention weights (B, N)."""
        h = self.encoder(X)                                   # (B, N, hidden)
        m = mask.unsqueeze(-1)                                # (B, N, 1)
        if self.pool == "attention":
            a = self.att_w(torch.tanh(self.att_V(h)) * torch.sigmoid(self.att_U(h)))  # (B,N,1)
            a = a.masked_fill(~m, NEG_INF)
            a = torch.softmax(a, dim=1)                       # over instances
            z = (a * h).sum(dim=1)                            # (B, hidden)
            aw = a.squeeze(-1)
        else:
            cnt = mask.sum(dim=1, keepdim=True).clamp(min=1)  # (B,1)
            z = (h * m).sum(dim=1) / cnt                      # (B, hidden)
            aw = mask.float() / cnt
        return self.head(z).squeeze(-1), aw


def pad_bags(bags):
    """list of (n_i, D) arrays -> (B, maxN, D) float tensor + (B, maxN) bool mask."""
    B = len(bags); D = bags[0].shape[1]
    maxN = max(b.shape[0] for b in bags)
    X = np.zeros((B, maxN, D), dtype=np.float32)
    mask = np.zeros((B, maxN), dtype=bool)
    for i, b in enumerate(bags):
        X[i, : b.shape[0]] = b
        mask[i, : b.shape[0]] = True
    return (torch.as_tensor(X, device=DEVICE),
            torch.as_tensor(mask, device=DEVICE))


def train_eval(bags_tr, y_tr, bags_va, y_va, in_dim, pool="attention",
               hidden=64, att_dim=32, dropout=0.25, lr=1e-3, wd=1e-3,
               max_epochs=300, patience=30, seed=0):
    torch.manual_seed(seed); np.random.seed(seed)
    model = MIL(in_dim, hidden, att_dim, pool, dropout).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=lr, weight_decay=wd)
    lossf = nn.BCEWithLogitsLoss()
    Xtr, mtr = pad_bags(bags_tr)
    ytr = torch.as_tensor(np.asarray(y_tr), dtype=torch.float32, device=DEVICE)
    Xva, mva = pad_bags(bags_va)
    yva = torch.as_tensor(np.asarray(y_va), dtype=torch.float32, device=DEVICE)

    best_val = np.inf; best_state = None; bad = 0
    for ep in range(max_epochs):
        model.train(); opt.zero_grad()
        logits, _ = model(Xtr, mtr)
        loss = lossf(logits, ytr)
        loss.backward(); opt.step()
        model.eval()
        with torch.no_grad():
            vloss = lossf(model(Xva, mva)[0], yva).item()
        if vloss < best_val - 1e-4:
            best_val = vloss; bad = 0
            best_state = {k: v.detach().clone() for k, v in model.state_dict().items()}
        else:
            bad += 1
            if bad >= patience:
                break
    if best_state is not None:
        model.load_state_dict(best_state)
    model.eval()
    with torch.no_grad():
        va_prob = torch.sigmoid(model(Xva, mva)[0]).cpu().numpy()
    return model, np.atleast_1d(va_prob)


def bag_attention(model, bag):
    """bag: (n, D) array -> attention weights (n,)."""
    model.eval()
    X, m = pad_bags([bag])
    with torch.no_grad():
        _, aw = model(X, m)
    return aw[0, : bag.shape[0]].cpu().numpy()
