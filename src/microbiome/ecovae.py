# Ecological Variational AutoEncoder

import torch
from torch import nn 
from torch.utils.data import Dataset, DataLoader 

import swanlab


class MicrobiomeDataset(Dataset):

    def __init__(self, abundance, sequences, augmentation: bool = True):
        self.abundance = abundance
        self.sequences = sequences
        self.augmentation = augmentation

    def __len__(self):
        return self.abundance.shape[0]

    def __getitem__(self, idx):
        abundance = self.abundance[idx]
        sequence = self.sequences[idx]

        if self.augmentation:
            abundance_aug, sequence_aug = self.taxa_permutation_invariance(
                abundance, sequence)
            return (abundance, sequence, abundance_aug, sequence_aug)

        return abundance, sequence

    @staticmethod
    def fake_datas(n_samples=1000, n_otus=256, seq_len=128):
        abundance = torch.rand(n_samples, n_otus)
        sequences = torch.randint(0, 4, (n_samples, n_otus, seq_len))

        return (abundance, sequences)

    @staticmethod
    def taxa_permutation_invariance(abundance, sequence):
        n_otus = abundance.shape[0]
        perm = torch.randperm(n_otus)
        abundance_aug = abundance[perm]
        sequence_aug = sequence[perm]

        return (abundance_aug, sequence_aug)


