

from typing import Literal, Any
from pydantic import BaseModel, Field
import numpy as np
import pandas as pd 
from matplotlib import pyplot as plt 
import seaborn as sns

from abc import ABC, abstractmethod

import torch
from torch import nn 
from torch.nn import functional as F
from torch.utils.data import Dataset, DataLoader 

import swanlab


# --- Dataset ---
## --- Preprocesse ---
class Preprocessor:
    '''Preprocessor for microbiome data'''

    def load_data(self, zOTU_fasta: str, zOTU_file: str, sintax_file: str, metainfo_file: str):
        """Read abundance and sequence data"""

        from Bio import SeqIO
        sequences: dict = {record.id: str(record.seq) for record in SeqIO.parse(zOTU_fasta, "fasta")}

        from microbiome import amplicon
        amplicon_operator = amplicon.Amplicon(
            otu_file_path=zOTU_file,
            sintax_file_path=sintax_file,
            metadata_file_path=metainfo_file
        )
        otu_table, taxonomy_table, metadata_table = amplicon_operator.features_parser()

        return (sequences,otu_table, taxonomy_table, metadata_table)
    
    def kmer(self, sequence, k):
        '''计算k-mer频率向量'''
        from itertools import product 
        bases = ['A', 'C', 'G', 'T']

    def _encode_sequences_with_hyenadna(self, sequences: pd.DataFrame, device: str = 'cpu', max_length: int = 256) -> torch.Tensor:
        from transformers import AutoModelForSequenceClassification, AutoTokenizer
        model_name = 'LongSafari/hyenadna-large-1m-seqlen'
        from transformers import AutoModel 
        model = AutoModel.from_pretrained(model_name, trust_remote_code=True)
        tokenizer = AutoTokenizer.from_pretrained(model_name, trust_remote_code=True)
        model.eval()
        model.to(device)
        tokens = tokenizer(sequences, return_tensors="pt", padding="max_length", truncation=True, max_length=max_length)
        return tokens # input_ids, attention_mask, etc.
    
    def _encode_sequences_with_nucleotide_transformer(self, sequences: pd.DataFrame, device: str = 'cpu', max_length: int = 256) -> torch.Tensor:
        pass 

    def _dna_bert(self):
        from transformers import AutoTokenizer, AutoModel

        tokenizer = AutoTokenizer.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)
        model = AutoModel.from_pretrained("zhihan1996/DNABERT-2-117M", trust_remote_code=True)


    def encode_sequences(self, sequences, method: Literal['demo', 'hyenadna', 'nucleotide_transformer'] = 'demo'):
        if method == 'hyenadna':
            return self._encode_sequences_with_hyenadna(sequences)
        elif method == 'nucleotide_transformer':
            return self._encode_sequences_with_nucleotide_transformer(sequences)
        else:
            raise ValueError("Unsupported encoding method")


## --- Data Augmentation ---
class Compose:
    '''组合多个数据增强或预处理操作'''

    def __init__(self, transforms):
        self.transforms = transforms

    def __call__(self, abundance: torch.Tensor, seq_features: torch.Tensor | None = None):
        if seq_features is not None:
            for transform in self.transforms:
                abundance, seq_features = transform(abundance, seq_features)
            return abundance, seq_features
        else:
            for transform in self.transforms:
                abundance = transform(abundance, seq_features=None)
            return abundance


class Permutation:
    '''顺序不变性增强：随机排列 OTU 顺序'''

    def __call__(self, abundance: torch.Tensor, seq_features: torch.Tensor | None = None):
        """随机排列 OTU 顺序"""
        n_otus: int = abundance.shape[0]
        perm: int = torch.randperm(n_otus)
        abundance_aug: torch.Tensor = abundance[perm]
        if seq_features is not None:
            seq_features_aug: torch.Tensor = seq_features[perm]
            return (abundance_aug, seq_features_aug)
        else:
            return abundance_aug


class NoiseInjection:
    '''噪声注入增强：在丰度向量中添加高斯噪声'''

    def __init__(self, noise_level: float = 0.01):
        self.noise_level = noise_level

    def __call__(self, abundance: torch.Tensor, seq_features: torch.Tensor | None = None):
        """在丰度向量中添加高斯噪声"""
        noise = torch.randn_like(abundance) * self.noise_level
        abundance_aug = abundance + noise
        abundance_aug = torch.clamp(abundance_aug, min=0.0)  # 保证非负
        if seq_features is not None:
            return (abundance_aug, seq_features)
        else:
            return abundance_aug


## --- Dataset ---
class MicrobiomeDataset(Dataset):
    '''Abundance + Sequence'''

    def __init__(self, abundance, seq_features=None):
        self.abundance = abundance              # (n_samples, n_otus)
        self.seq_features = seq_features        # (n_otus, feat_dim) or None

    def __len__(self):
        return self.abundance.shape[0]

    def __getitem__(self, idx):
        abund = self.abundance[idx]
        if self.seq_features is not None:
            seq_feat = self.seq_features  # 返回整个特征矩阵，或 idx 对应的？此处应该返回所有OTU的特征（因为重建需要所有OTU）
            # 通常 seq_features 是 (n_otus, d)，每个样本共享，所以直接返回即可
            return abund, seq_feat
        else:
            return abund
        
    def permute_features(self, perm):
        """按 perm 重排丰度矩阵的列（特征维度）"""

        self.abundance = self.abundance[:, perm]
        # 如果 seq_features 也存在且需要同步置换，按同样方式处理
        if self.seq_features is not None:
            self.seq_features = self.seq_features[:, perm]
        

class MicrobiomeDatasetSampleAugmentation(Dataset):
    """微生物组数据集： Abundance + Sequence Features
    对每个__getitem__的数据都进行乱序。
    """

    def __init__(self, abundance, seq_features = None, transform = None):
        """
        abundance: (n_samples, n_otus)  相对丰度（或转换后的丰度）
        seq_features: (n_otus, feature_dim)  每个OTU的预提取嵌入（如DNABERT/HyenaDNA的输出）
        transform: 数据转换函数
        """
        super().__init__()
        self.abundance = abundance
        self.seq_features = seq_features
        self.transform = transform

    def __len__(self):
        return self.abundance.shape[0]

    def __getitem__(self, idx):
        abund = self.abundance[idx]            # [n_otus]
        if self.transform is not None:
            if self.seq_features is not None:
                # with sequence
                abund_aug, seq_features_aug = self.transform(abund, self.seq_features)
                return (abund, self.seq_features, abund_aug, seq_features_aug)
            else:
                # without sequence
                abund_aug = self.transform(abund, None)
                return (abund, abund_aug)
        else:
            if self.seq_features is not None:
                # 如果没有transform，直接返回原始数据
                return (abund, self.seq_features)
            else:
                return abund

    @staticmethod
    def fake_datas(n_samples=1000, n_otus=5000, feature_dim=768) -> tuple[torch.Tensor, torch.Tensor]:
        """生成虚拟数据用于测试"""
        abundance = torch.rand(n_samples, n_otus)
        abundance = abundance / abundance.sum(dim=1, keepdim=True)   # 归一化为相对丰度
        seq_features = torch.randn(n_otus, feature_dim)              # 假装是预训练嵌入
        return (abundance, seq_features)


# --- Model ---
class DeepAVAE(nn.Module):
    '''Abundance VAE: 仅使用丰度向量进行编码和解码'''

    def __init__(self, input_dim, latent_dim):
        super().__init__()
        # 编码器：加隐藏层 + BatchNorm + 激活
        hiddent_dim: int = 128
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, hiddent_dim),
            nn.BatchNorm1d(hiddent_dim),
            nn.ReLU(),
            nn.Linear(hiddent_dim, hiddent_dim*2),
            nn.BatchNorm1d(hiddent_dim*2),
            nn.ReLU()
        )
        self.mu_head = nn.Linear(hiddent_dim*2, latent_dim)
        self.logvar_head = nn.Linear(hiddent_dim*2, latent_dim)

        # 解码器：从 latent 到 hidden，再逐步放大
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, hiddent_dim*2),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(hiddent_dim*2, hiddent_dim),
            nn.ReLU(),
            nn.Linear(hiddent_dim, input_dim),
            # nn.Softmax(dim=-1)  # 如果你确定目标已归一化
        )

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        h = self.encoder(x)
        mu = self.mu_head(h)
        logvar = self.logvar_head(h)
        z = self.reparameterize(mu, logvar)
        logits = self.decoder(z)
        return (logits, mu, logvar)
    

class DeepASVAE(nn.Module):
    '''Abundance + Sequence VAE: 使用丰度向量和序列嵌入进行编码和解码'''

    def __init__(self, input_dim, latent_dim, n_otu):
        super().__init__()
        # 编码器：加隐藏层 + BatchNorm + 激活
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.ReLU()
        )
        self.mu_head = nn.Linear(256, latent_dim)
        self.logvar_head = nn.Linear(256, latent_dim)

        # 解码器：从 latent 到 hidden，再逐步放大
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 256),
            nn.ReLU(),
            nn.Dropout(0.2),
            nn.Linear(256, 512),
            nn.ReLU(),
            nn.Linear(512, n_otu),
            # nn.Softmax(dim=-1)  # 如果你确定目标已归一化
        )

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x):
        h = self.encoder(x)
        mu = self.mu_head(h)
        logvar = self.logvar_head(h)
        z = self.reparameterize(mu, logvar)
        logits = self.decoder(z)
        return (logits, mu, logvar)
    

class VAMB_VAE(nn.Module):
    """
    仿 VAMB 架构的多模态 VAE：
    输入 = 丰度向量 + 序列特征向量（拼接）
    编码器和解码器均使用两层 MLP，潜在维度 32~64
    """

    def __init__(self, abund_dim, seq_dim, latent_dim=32, beta=1.0):
        super().__init__()
        self.input_dim = abund_dim + seq_dim
        self.latent_dim = latent_dim
        self.beta = beta

        # 编码器：512 → 256（与 VAMB 一致）
        self.encoder = nn.Sequential(
            nn.Linear(self.input_dim, 512),
            nn.BatchNorm1d(512),
            nn.ReLU(),
            nn.Linear(512, 256),
            nn.BatchNorm1d(256),
            nn.ReLU()
        )
        self.mu_head = nn.Linear(256, latent_dim)
        self.logvar_head = nn.Linear(256, latent_dim)

        # 解码器：256 → 512 → input_dim（对称）
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 512),
            nn.ReLU(),
            nn.Linear(512, self.input_dim)
        )

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def forward(self, x_abund, x_seq):
        # 拼接丰度与序列特征
        x = torch.cat([x_abund, x_seq], dim=1)
        h = self.encoder(x)
        mu = self.mu_head(h)
        logvar = self.logvar_head(h)
        z = self.reparameterize(mu, logvar)
        recon = self.decoder(z)  # 重建整个拼接向量
        # 拆分开重建的丰度和序列部分
        recon_abund = recon[:, :x_abund.shape[1]]
        recon_seq = recon[:, x_abund.shape[1]:]
        return recon_abund, recon_seq, mu, logvar

    def loss(self, x_abund, x_seq, recon_abund, recon_seq, mu, logvar):
        # 重建损失：分别计算丰度和序列特征的 MSE，可赋予不同权重
        recon_loss_abund = F.mse_loss(recon_abund, x_abund)
        recon_loss_seq = F.mse_loss(recon_seq, x_seq)
        recon_loss = recon_loss_abund + recon_loss_seq
        # KL 散度
        kl = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1).mean()
        return recon_loss + self.beta * kl
    

# --- Training ---
class Trainer(ABC):
    '''Abstract base class for training models'''

    @abstractmethod
    def train(self, model, dataloader, optimizer, scheduler, device, epochs=10, beta=0.1, **kwargs): ...

    @abstractmethod
    def inference(self, **kswargs): ...

    def cluster_metric(self, latent_representations, cluster_range=[1, 10]):
        from sklearn.metrics import silhouette_score
        from sklearn.cluster import KMeans

        if isinstance(latent_representations, torch.Tensor):
            latent_representations = latent_representations.cpu()
        
        cluster_num: list = []
        scopes: list = []
        for i in range(cluster_range[0], cluster_range[1]):
            # clustering
            model = KMeans(n_clusters=i, random_state=42)
            labels = model.fit_predict(latent_representations)
            # if len(labels) < 2:
            #     cluster_num.append(i)
            #     scopes.append(0)
            # eval
            scope = silhouette_score(latent_representations, labels)
            cluster_num.append(i)
            scopes.append(scope)

        df = pd.DataFrame(
            {
                'cluster': cluster_num,
                'scope': scopes,
            }
        )

        # 最大聚集系数对应的聚类书
        max_scope = df.sort_values(by='scope', ascending=False).iloc[0]['cluster']

        return max_scope.astype(int)

    def save_model(
        self,
        model: nn.Module,
        path: str,
        method: Literal["state_dict", "jit"] = "state_dict",
        jit_model: Literal['script', 'trace'] = 'script',               # True: script, False: trace
        example_input: torch.Tensor = None,   # for trace
    ):
        '''Save final trained model'''

        # 先剥离分布式wrapper，获取原始模块
        raw_model = model
        if isinstance(model, (nn.DataParallel, nn.parallel.DistributedDataParallel)):
            raw_model = model.module

        if method == "state_dict":
            torch.save(raw_model.state_dict(), path)
        elif method == "jit":
            # 转为eval模式，避免训练行为固化
            raw_model.eval()
            if jit_model == 'script':
                scripted_model = torch.jit.script(raw_model)
            elif jit_model == 'trace':
                if example_input is None:
                    raise ValueError("example_input must be provided when jit_model='trace'")
                scripted_model = torch.jit.trace(raw_model, example_input)
            else:
                raise ValueError(f"Unsupported use_script value: {use_script}. Use 'script' or 'trace'.")
            torch.jit.save(scripted_model, path)
        else:
            raise ValueError(f"Unsupported saving method: {method}. Use 'state_dict' or 'jit'.")

        print(f"Model saved to {path} using method '{method}'.")

    def load_model(
            self,
            path: str,
            method: Literal["state_dict", "jit"] = "jit",
            device: str = 'cpu'
    ):
        '''Load trained model'''
        if method == "state_dict":
            state_dict = torch.load(path, map_location=torch.device(device))
            return state_dict 
        
        elif method == "jit":
            model = torch.jit.load(path, map_location=torch.device(device))
            return model
        
        else:
            raise ValueError(f"Unsupported loading method: {method}. Use 'state_dict' or 'jit'.")
    

class ATrainer(Trainer):
    '''Trainer for Abundance'''

    def loss_fn(self, x, logits, mu, logvar, beta): 
        # 重构损失：交叉熵
        # recon:解码器最后一层没有经过 Softmax 的输出的 logits
        log_prob = F.log_softmax(logits, dim=-1)          # 稳定 log 概率
        recon_loss = -torch.sum(x * log_prob, dim=-1).mean()

        kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=-1).mean()
        loss = recon_loss + beta * kl_loss

        return (loss, recon_loss, kl_loss)

    def train(self, model, dataloader, test_loader, optimizer, scheduler, device, epochs, beta, **kwargs) -> None:
        for epoch in range(epochs):
            model.train()
            total_loss = 0.0
            total_recon_loss = 0.0
            total_kl_loss = 0.0

            for batch in dataloader:
                # --- without augmentation ---
                abund = batch
                abund = abund.to(device)
                recon, mu, logvar = model(abund)
                loss, recon_loss, kl_loss = self.loss_fn(
                    abund, recon, mu, logvar, beta
                )

                # --- with augmentation ---
                # abund, abund_aug = batch
                # abund = abund.to(device)
                # abund_aug = abund_aug.to(device)
                # recon, mu, logvar = model(abund_aug)
                # loss, recon_loss, kl_loss = self.loss_fn(
                #     abund_aug, recon, mu, logvar, beta
                # )

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                total_loss += loss.item()
                total_recon_loss += recon_loss.item()
                total_kl_loss += kl_loss.item()

            if (epoch + 1) % 5 == 0:
                batch_num: int = len(dataloader)
                metrics: dict = {
                    'epoch': epoch + 1,
                    'current_lr_opti': optimizer.param_groups[0]['lr'],
                    't_loss': total_loss / batch_num,
                    't_recon': total_recon_loss / batch_num,
                    't_kl': total_kl_loss / batch_num,
                    'max_scope': self.validate(model=model, test_loader=test_loader, device=device)
                }
                content: str = ''
                for key, value in metrics.items():
                    content += f'{key}: {value:.9f} | '
                print(content)

            scheduler.step()

        print("训练完成。")
        return model

    def inference(self, model, test_loader, device):
        """给定丰度向量和全局嵌入，返回样本的潜在表征"""

        model.eval()
        with torch.no_grad():
            mus = []
            for abund in test_loader:
                abund = abund.to(device)
                _, mu, _ = model(abund)
                mus.append(mu)
            return torch.cat(mus)  # 返回潜在表征
        
    def validate(self, model, test_loader, device):
        '''验证'''

        mu = self.inference(model, test_loader, device)
        max_scope = self.cluster_metric(latent_representations=mu, cluster_range=[2, 10])

        return max_scope


class ASTrainer(Trainer):
    '''Trainer for Abundance and Sequence Featrure'''

    def loss_fn(self, x, logits, mu, logvar, beta): 
        # 重构损失：交叉熵
        # recon:解码器最后一层没有经过 Softmax 的输出的 logits
        log_prob = F.log_softmax(logits, dim=-1)          # 稳定 log 概率
        recon_loss = -torch.sum(x * log_prob, dim=-1).mean()

        kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=-1).mean()
        loss = recon_loss + beta * kl_loss

        return (loss, recon_loss, kl_loss)

    def train(self, model, dataloader, test_loader, optimizer, scheduler, device, epochs=10, beta=0.1, **kwargs) -> None:
        for epoch in range(epochs):
            model.train()
            total_loss = 0.0
            total_recon_loss = 0.0
            total_kl_loss = 0.0

            for batch in dataloader:
                # batch 返回: (abund, seq_features_global, abund_aug, seq_features_aug)
                abund_orig, _, abund_aug, seq_aug = batch
                # 使用增强后的数据构建样本表征，并作为 VAE 输入
                abund_aug = abund_aug.to(device)
                seq_aug = seq_aug.to(device)
                sample_feat = torch.bmm(abund_aug.unsqueeze(1), seq_aug).squeeze(1)

                recon, mu, logvar = model(sample_feat)
                loss, recon_loss, kl_loss = self.loss_fn(abund_aug, recon, mu, logvar, beta)

                optimizer.zero_grad()
                loss.backward()
                optimizer.step()

                total_loss += loss.item()
                total_recon_loss += recon_loss.item()
                total_kl_loss += kl_loss.item()

            if (epoch + 1) % 5 == 0:
                # current_lr_sche = scheduler.get_last_lr()[0]
                batch_num: int = len(dataloader)
                metrics: dict = {
                    'epoch': epoch + 1,
                    'current_lr_opti': optimizer.param_groups[0]['lr'],
                    't_loss': total_loss / batch_num,
                    't_recon': total_recon_loss / batch_num,
                    't_kl': total_kl_loss / batch_num,
                    'max_scope': self.validate(model=model, test_loader=test_loader, device=device)
                }
                content: str = ''
                for key, value in metrics.items():
                    content += f'{key}: {value:.9f} | '
                print(content)
            
            scheduler.step()

        print("训练完成。")
        return model

    def inference(self, model, test_loader, device):
        """给定丰度向量和全局嵌入，返回样本的潜在表征"""

        model.eval()
        with torch.no_grad():
            mus = []
            for abund, seq in test_loader:
                abund = abund.to(device)
                seq = seq.to(device)
                sample_feat = torch.bmm(abund.unsqueeze(1), seq).squeeze(1)
                # h = model.encoder(sample_feat)
                _, mu, _ = model(sample_feat)
                mus.append(mu)
            return torch.cat(mus)  # 返回潜在表征
        
    def validate(self, model, test_loader, device):
        '''验证'''

        mu = self.inference(model, test_loader, device)
        max_scope = self.cluster_metric(latent_representations=mu, cluster_range=[2, 10])

        return max_scope
    

class ASTrainerEpochAugmentation(Trainer):

    # def loss_fn(self, x, logits, mu, logvar, beta): 
    def loss_fn(self, x, logits, mu, logvar, beta): 
        # 重构损失：交叉熵
        # recon:解码器最后一层没有经过 Softmax 的输出的 logits
        log_prob = F.log_softmax(logits, dim=-1)          # 稳定 log 概率
        recon_loss = -torch.sum(x * log_prob, dim=-1).mean()

        kl_loss = -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=-1).mean()
        loss = recon_loss + beta * kl_loss

        return (loss, recon_loss, kl_loss)

    def train(
            self, abundance, seq_features, device='cuda:0', epochs=10, beta=0.1, 
            augmentation: bool = False,
            **kwargs
        ) -> None:
        
        # --- repeat ---
        torch.manual_seed(42)
        torch.cuda.manual_seed(42)
        torch.cuda.manual_seed_all(42)

        # --- split dataset ---
        from sklearn.model_selection import train_test_split
        train_abundance, test_abundance = train_test_split(
            abundance, test_size=0.2, random_state=123,
        )

        loader_params = {
            'batch_size': 128,
            'num_workers': 8,
            'persistent_workers': True,
            'pin_memory': True,
        }
        train_dataset = MicrobiomeDataset(abundance=train_abundance, seq_features=seq_features)
        train_loader = DataLoader(train_dataset, shuffle=True, **loader_params)
        test_dataset = MicrobiomeDataset(abundance=test_abundance, seq_features=seq_features)
        test_loader = DataLoader(test_dataset, shuffle=False, **loader_params)

        input_dim = abundance.shape[1] # n_otus
        model = DeepAVAE(input_dim=input_dim, latent_dim=32).to(device)
        # model = DeepASVAE().to(device)

        optimizer = torch.optim.AdamW(params=model.parameters(), lr=1e-2)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

        for epoch in range(epochs):
            if augmentation:
                # Augmentation per EPOCH. 
                perm = torch.randperm(input_dim)
                train_loader.dataset.permute_features(perm)
                test_loader.dataset.permute_features(perm)
                # print(f'Augmentation: [{epoch + 1}]')

            model.train()
            total_loss = 0.0
            total_recon_loss = 0.0
            total_kl_loss = 0.0

            for batch in train_loader:
                abund = batch.to(device)
                logits, mu, logvar = model(abund)
                loss, recon_loss, kl_loss = self.loss_fn(abund, logits, mu, logvar, beta)

                optimizer.zero_grad()
                loss.backward()
                total_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()

                total_loss += loss.item()
                total_recon_loss += recon_loss.item()
                total_kl_loss += kl_loss.item()

            if (epoch + 1) % 50 == 0:
                batch_num: int = len(train_loader)
                metrics: dict = {
                    'epoch': epoch + 1,
                    'current_lr_opti': optimizer.param_groups[0]['lr'],
                    't_loss': total_loss / batch_num,
                    't_recon': total_recon_loss / batch_num,
                    't_kl': total_kl_loss / batch_num,
                    'max_scope': self.validate(model=model, test_loader=test_loader, device=device)
                }
                content: str = ''
                for key, value in metrics.items():
                    content += f'{key}: {value:.9f} | '
                print(content)
            
            if scheduler:
                scheduler.step()

        print("训练完成。")
        return model, test_loader

    def inference(self, model, test_loader, device):
        """给定丰度向量和全局嵌入，返回样本的潜在表征"""

        model.eval()
        with torch.no_grad():
            mus = []
            for batch in test_loader:
                # abund = abund.to(device)
                # seq = seq.to(device)
                # sample_feat = torch.bmm(abund.unsqueeze(1), seq).squeeze(1)
                # _, mu, _ = model(sample_feat)
                abund = batch.to(device)
                _, mu, _ = model(abund)
                mus.append(mu)
            return torch.cat(mus)  # 返回潜在表征
        
    def validate(self, model, test_loader, device):
        '''验证'''

        mu = self.inference(model, test_loader, device)
        max_scope = self.cluster_metric(latent_representations=mu, cluster_range=[2, 10])

        return max_scope
    

class ASTrainerEpochAugmentationVAMB(Trainer):

    def train(
            self, abundance, seq_feature, device='cuda:0', epochs=10, beta=0.1, 
            augmentation: bool = False,
            **kwargs
        ) -> None:
        
        # --- repeat ---
        torch.manual_seed(42)
        torch.cuda.manual_seed(42)
        torch.cuda.manual_seed_all(42)

        # --- split dataset ---
        from sklearn.model_selection import train_test_split
        train_abundance, test_abundance = train_test_split(
            abundance, test_size=0.2, random_state=123,
        )

        loader_params = {
            'batch_size': 128,
            'num_workers': 8,
            'persistent_workers': True,
            'pin_memory': True,
        }
        train_dataset = MicrobiomeDataset(abundance=train_abundance, seq_features=seq_feature)
        train_loader = DataLoader(train_dataset, shuffle=True, **loader_params)
        test_dataset = MicrobiomeDataset(abundance=test_abundance, seq_features=seq_feature)
        test_loader = DataLoader(test_dataset, shuffle=False, **loader_params)

        input_dim = abundance.shape[1] # n_otus
        model = DeepAVAE(input_dim=input_dim, latent_dim=32).to(device)
        # model = DeepASVAE().to(device)

        optimizer = torch.optim.AdamW(params=model.parameters(), lr=1e-2)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)

        for epoch in range(epochs):
            if augmentation:
                # Augmentation per EPOCH. 
                perm = torch.randperm(input_dim)
                train_loader.dataset.permute_features(perm)
                test_loader.dataset.permute_features(perm)
                # print(f'Augmentation: [{epoch + 1}]')

            model.train()
            total_loss = 0.0
            total_recon_loss = 0.0
            total_kl_loss = 0.0

            for batch in train_loader:
                abund = batch.to(device)
                logits, mu, logvar = model(abund)
                loss, recon_loss, kl_loss = self.loss_fn(abund, logits, mu, logvar, beta)

                optimizer.zero_grad()
                loss.backward()
                total_norm = torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=1.0)
                optimizer.step()

                total_loss += loss.item()
                total_recon_loss += recon_loss.item()
                total_kl_loss += kl_loss.item()

            if (epoch + 1) % 50 == 0:
                batch_num: int = len(train_loader)
                metrics: dict = {
                    'epoch': epoch + 1,
                    'current_lr_opti': optimizer.param_groups[0]['lr'],
                    't_loss': total_loss / batch_num,
                    't_recon': total_recon_loss / batch_num,
                    't_kl': total_kl_loss / batch_num,
                    'max_scope': self.validate(model=model, test_loader=test_loader, device=device)
                }
                content: str = ''
                for key, value in metrics.items():
                    content += f'{key}: {value:.9f} | '
                print(content)
            
            if scheduler:
                scheduler.step()

        print("训练完成。")
        return model, test_loader

    def inference(self, model, test_loader, device):
        """给定丰度向量和全局嵌入，返回样本的潜在表征"""

        model.eval()
        with torch.no_grad():
            mus = []
            for batch in test_loader:
                # abund = abund.to(device)
                # seq = seq.to(device)
                # sample_feat = torch.bmm(abund.unsqueeze(1), seq).squeeze(1)
                # _, mu, _ = model(sample_feat)
                abund = batch.to(device)
                _, mu, _ = model(abund)
                mus.append(mu)
            return torch.cat(mus)  # 返回潜在表征
        
    def validate(self, model, test_loader, device):
        '''验证'''

        mu = self.inference(model, test_loader, device)
        max_scope = self.cluster_metric(latent_representations=mu, cluster_range=[2, 10])

        return max_scope


# --- Analysis ---
class Analyzer:
    '''The downstream analysis module for latent representations'''

    def reducer(
            self, 
            latent_representations, 
            method: Literal['umap', 'tsne', 'pca'] = 'umap', 
            n_components: int = 2,
    ):
        """对潜在表征进行降维"""

        if method == 'umap':
            import umap
            reducer = umap.UMAP(n_components=n_components, random_state=42)
            reduced = reducer.fit_transform(latent_representations.cpu().numpy())
            return reduced
        elif method == 'tsne':
            from sklearn.manifold import TSNE
            tsne = TSNE(n_components=n_components, random_state=42)
            reduced = tsne.fit_transform(latent_representations.cpu().numpy())
            return reduced
        elif method == 'pca':
            from sklearn.decomposition import PCA
            pca = PCA(n_components=n_components)
            reduced = pca.fit_transform(latent_representations.cpu().numpy())
            return reduced
        else:
            raise ValueError("Unsupported dimensionality reduction method")

    def cluster(
        self,
        latent_representations,          # torch.Tensor 或 np.ndarray
        method: Literal['kmeans', 'hierarchical', 'dbscan', 'optics',
                        'gmm', 'spectral', 'affinity_propagation',
                        'mean_shift', 'birch'] = 'kmeans',
        n_clusters: int = 5,            # 部分算法会忽略此参数
        **kwargs: Any                   # 传递算法专属参数（如 eps, min_samples 等）
    ) -> np.ndarray:
        """
        对潜在表征进行聚类，支持多种算法。

        参数
        ----------
        latent_representations : torch.Tensor 或 np.ndarray
            形状为 (n_samples, n_features) 的潜在表示。
        method : str
            聚类算法标识。
        n_clusters : int, 默认 5
            簇数（仅对需要预设簇数的算法有效）。
        **kwargs : dict
            传递给具体聚类器的其他参数，例如：
            - DBSCAN/OPTICS: eps, min_samples
            - MeanShift: bandwidth, bin_seeding
            - Spectral: affinity, gamma, assign_labels
            - GMM: covariance_type, init_params, reg_covar
            - Birch: threshold, branching_factor

        返回
        -------
        labels : np.ndarray, shape (n_samples,)
            聚类标签，噪声点标记为 -1（如果算法支持）。
        """
        # 统一转为 numpy 数组（兼容 PyTorch Tensor）
        from sklearn.cluster import (KMeans, AgglomerativeClustering, DBSCAN,
                                    OPTICS, SpectralClustering, AffinityPropagation,
                                    MeanShift, estimate_bandwidth, Birch)
        from sklearn.mixture import GaussianMixture

        if hasattr(latent_representations, 'cpu'):  # torch.Tensor
            X = latent_representations.detach().cpu().numpy()
        else:
            X = np.asarray(latent_representations)

        # 根据 method 选择并实例化聚类器
        if method == 'kmeans':
            model = KMeans(n_clusters=n_clusters, random_state=42, **kwargs)
            labels = model.fit_predict(X)

        elif method == 'hierarchical':
            # AgglomerativeClustering 需要 n_clusters 或 distance_threshold
            # 这里默认使用 n_clusters，若用户传入 distance_threshold 则优先
            params = {'n_clusters': n_clusters}
            if 'distance_threshold' in kwargs:
                params['n_clusters'] = None          # n_clusters 与 distance_threshold 互斥
                params['distance_threshold'] = kwargs.pop('distance_threshold')
            model = AgglomerativeClustering(**params, **kwargs)
            labels = model.fit_predict(X)

        elif method == 'dbscan':
            # 默认参数：eps=0.5, min_samples=5，可通过 kwargs 覆盖
            model = DBSCAN(eps=kwargs.pop('eps', 0.5),
                           min_samples=kwargs.pop('min_samples', 5), **kwargs)
            labels = model.fit_predict(X)

        elif method == 'optics':
            # OPTICS 默认 min_samples=5，可通过 kwargs 调整
            model = OPTICS(min_samples=kwargs.pop('min_samples', 5), **kwargs)
            labels = model.fit_predict(X)

        elif method == 'gmm':
            # GaussianMixture 返回概率，取 argmax 得到硬标签
            model = GaussianMixture(n_components=n_clusters,
                                    random_state=42, **kwargs)
            model.fit(X)
            labels = model.predict(X)        # 也可用 model.fit_predict(X)

        elif method == 'spectral':
            # 谱聚类，默认使用 RBF 核，n_clusters 必须提供
            model = SpectralClustering(n_clusters=n_clusters,
                                       random_state=42,
                                       assign_labels='kmeans', **kwargs)
            labels = model.fit_predict(X)

        elif method == 'affinity_propagation':
            # 自动决定簇数，忽略 n_clusters
            # 关键参数：damping, preference，可通过 kwargs 传入
            model = AffinityPropagation(random_state=42, **kwargs)
            labels = model.fit_predict(X)

        elif method == 'mean_shift':
            # 自动决定簇数，忽略 n_clusters
            # 带宽可通过 kwargs 传入或自动估计
            bandwidth = kwargs.pop('bandwidth', None)
            if bandwidth is None:
                bandwidth = estimate_bandwidth(X, n_jobs=-1)
            model = MeanShift(bandwidth=bandwidth, bin_seeding=True, **kwargs)
            labels = model.fit_predict(X)

        elif method == 'birch':
            # Birch 需要 n_clusters（也可结合全局聚类器）
            model = Birch(n_clusters=n_clusters, **kwargs)
            labels = model.fit_predict(X)

        else:
            raise ValueError(f"不支持的聚类方法: {method}")

        return labels

    def clustering_with_scope(self, latent_representations, method: str = 'kmeans', cluster_range: list = [2, 10]):
        from sklearn.metrics import silhouette_score

        if isinstance(latent_representations, torch.Tensor):
            latent_representations = latent_representations.cpu()
        
        cluster_num: list = []
        scopes: list = []
        for i in range(cluster_range[0], cluster_range[1]):
            labels = self.cluster(latent_representations, method, n_clusters=i)
            scope = silhouette_score(latent_representations, labels)
            cluster_num.append(i)
            scopes.append(scope)

        df = pd.DataFrame(
            {
                'cluster': cluster_num,
                'scope': scopes,
            }
        )

        # 最大聚集系数对应的聚类书
        max_scope = df.sort_values(by='scope', ascending=False).iloc[0]['cluster']

        return (df, max_scope.astype(int))
    

def clustering_curve(mu, method):
    analyzer = Analyzer()
    
    # --- reduce ---
    fig = plt.figure(figsize=(3, 3))
    df, max_scope = analyzer.clustering_with_scope(latent_representations=mu, method=method)
    sns.barplot(data=df, x='cluster', y='scope')    

    return max_scope


def clustering(mu, cluster_num):
    analyzer = Analyzer()
    
    # --- reduce ---
    reduced_mu_umap = analyzer.reducer(latent_representations=mu, method='umap', n_components=2)
    reduced_mu_tsne = analyzer.reducer(latent_representations=mu, method='tsne', n_components=2)
    reduced_mu_pca = analyzer.reducer(latent_representations=mu, method='pca', n_components=2)
    labels_kmeans = analyzer.cluster(latent_representations=mu, method='kmeans', n_clusters=cluster_num)
    # labels_dbscan = analyzer.cluster(latent_representations=mu, method='dbscan', eps=0.5, min_samples=5)
    
    # --- dataframe ---
    reduced_mu_umap_df = pd.DataFrame(reduced_mu_umap, columns=['UMAP1', 'UMAP2'])
    reduced_mu_tsne_df = pd.DataFrame(reduced_mu_tsne, columns=['tSNE1', 'tSNE2'])  
    reduced_mu_pca_df = pd.DataFrame(reduced_mu_pca, columns=['PCA1', 'PCA2'])
    reduced_mu_umap_df['kmeans'] = labels_kmeans
    reduced_mu_tsne_df['kmeans'] = labels_kmeans
    reduced_mu_pca_df['kmeans'] = labels_kmeans
    # reduced_mu_umap_df['dbscan'] = labels_dbscan
    # reduced_mu_tsne_df['dbscan'] = labels_dbscan
    # reduced_mu_pca_df['dbscan'] = labels_dbscan

    # --- visualization ---
    cols: int = 3
    rows: int = 1
    fig, axs = plt.subplots(rows, cols, figsize=(cols * 6, rows * 4))
    i = 0
    sns.scatterplot(data=reduced_mu_umap_df, x='UMAP1', y='UMAP2', hue='kmeans', palette='tab10', ax=axs[i])
    axs[i].set_title('UMAP of Latent Representations')
    i = 1
    sns.scatterplot(data=reduced_mu_tsne_df, x='tSNE1', y='tSNE2', hue='kmeans', palette='tab10', ax=axs[i])
    axs[i].set_title('tSNE of Latent Representations')
    i = 2
    sns.scatterplot(data=reduced_mu_pca_df, x='PCA1', y='PCA2', hue='kmeans', palette='tab10', ax=axs[i])
    axs[i].set_title('PCA of Latent Representations')

    return None


class Footer:
    def __init__(self):
        pass