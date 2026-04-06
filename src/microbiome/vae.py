
from loguru import logger
from pathlib import Path

import torch
import torch.nn as nn
# import torchmetrics

import numpy as np
import pandas as pd 

import sklearn
import skbio


class EncoderMLP(nn.Module):
    def __init__(self, input_dim, latent_dim, dropout: float = 0.1):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Linear(input_dim, 4096),
            nn.BatchNorm1d(4096),
            nn.ReLU(),
            nn.Dropout(dropout) if dropout > 0 else nn.Identity(),
            nn.Linear(4096, 1024),
            nn.ReLU(),
            nn.Linear(1024, 256),
            nn.ReLU(),
        )
        self.mu = nn.Linear(256, latent_dim)
        self.logvar = nn.Linear(256, latent_dim)

    def forward(self, x):
        h = self.encoder(x)
        mu = self.mu(h)
        logvar = self.logvar(h)

        return mu, logvar


class DecoderMLP(nn.Module):
    def __init__(self, latent_dim, output_dim):
        super().__init__()
        self.decoder = nn.Sequential(
            nn.Linear(latent_dim, 256),
            nn.ReLU(),
            nn.Linear(256, 1024),
            nn.ReLU(),
            nn.Linear(1024, 4096),
            nn.ReLU(),
            nn.Linear(4096, output_dim),
        )

    def forward(self, z):
        y_hat = self.decoder(z)

        return y_hat


class VAE(nn.Module):
    '''Variational Autoencoder for dimensionality reduction of microbiome data.'''

    def __init__(self, input_dim, latent_dim=8):
        super().__init__()
        self.encoder = EncoderMLP(input_dim, latent_dim)
        self.decoder = DecoderMLP(latent_dim, input_dim)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std) # eps ~ N(0, I)
        return mu + eps * std

    def forward(self, x):
        mu, logvar = self.encoder(x)
        z = self.reparameterize(mu, logvar)
        y_hat = self.decoder(z)
        
        return y_hat, mu, logvar


# class VAE(nn.Module):
#     '''Variational Autoencoder for dimensionality reduction of microbiome data.'''

#     def __init__(self, input_dim, latent_dim=8):
#         super().__init__()

#         self.encoder = nn.Sequential(
#             nn.Linear(input_dim, 256),
#             nn.ReLU(),
#             nn.Linear(256, 64),
#             nn.ReLU()
#         )
#         self.mu = nn.Linear(64, latent_dim)
#         self.logvar = nn.Linear(64, latent_dim)

#         self.decoder = nn.Sequential(
#             nn.Linear(latent_dim, 64),
#             nn.ReLU(),
#             nn.Linear(64, 256),
#             nn.ReLU(),
#             nn.Linear(256, input_dim)
#         )

#     def reparameterize(self, mu, logvar):
#         std = torch.exp(0.5 * logvar)
#         eps = torch.randn_like(std) # eps ~ N(0, I)

#         return mu + eps * std

#     def forward(self, x):
#         h = self.encoder(x)
#         mu = self.mu(h)
#         logvar = self.logvar(h)
#         z = self.reparameterize(mu, logvar)

#         return self.decoder(z), mu, logvar


class VAEPipeline:
    '''Pipeline for training a VAE on microbiome OTU data.'''
    
    def __init__(self, seed: int = 123):
        logger.info(f'Reproducibility: Setting random seed to {seed} for torch, numpy, and random.')
        torch.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
        # np.random.seed(seed)
        # import random
        # random.seed(seed)

    @staticmethod
    def preprocess(otu_table: pd.DataFrame) -> pd.DataFrame:
        '''Preprocess the OTU table by adding a pseudocount, applying CLR transformation, and scaling the data.
        Args:
            otu_table: OTU table with OTU as rows and samples as columns.
        Returns:
            Preprocessed OTU table with samples as rows and OTUs as columns.
        '''

        logger.info('Adding pseudocount and CLR transformation to OTU table...')
        clr_out_table: pd.DataFrame = pd.DataFrame(
            data= skbio.stats.composition.clr(
                skbio.stats.composition.multi_replace(otu_table.T), # sample X OTU_ID, with pseudocount that replaces zeros.
            ),
            index= otu_table.T.index,
            columns= otu_table.T.columns,
        )

        logger.info('Scale [SKLEARN StandardScaler] CLR-transformed OTU table...')
        scaler = sklearn.preprocessing.StandardScaler()
        clr_out_table[clr_out_table.columns] = scaler.fit_transform(clr_out_table[clr_out_table.columns])

        return clr_out_table

    def load_dataset(self, otu_table: pd.DataFrame, batch_size: int = 32, num_workers: int = 4):
        '''Load the OTU table as a PyTorch dataset.'''
        
        X = torch.tensor(otu_table.values, dtype=torch.float32)
        dataset = torch.utils.data.TensorDataset(X)
        loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True, pin_memory=True, num_workers=num_workers)

        return loader

    def loss_function(self, recon_x, x, mu, logvar, beta: float = 1.0):
        '''Calculate the VAE loss, which is a combination of reconstruction loss and KL divergence.'''

        recon_loss = torch.nn.MSELoss()(recon_x, x)
        kl = -0.5 * torch.mean(1 + logvar - mu.pow(2) - logvar.exp())
        return recon_loss + beta * kl, recon_loss, kl

    def fit(self, model_configs: dict, otu_table: pd.DataFrame, 
            device: str = 'cpu', 
            epochs: int = 10, batch_size: int = 32, lr: float = 1e-3, 
            ckpt_dir: str = './logs'
        ):
        '''Fit the VAE model to the OTU table.'''

        ckpt_dir: Path = Path(ckpt_dir)
        if not ckpt_dir.exists():
            ckpt_dir.mkdir(parents=True, exist_ok=True)   
            logger.info(f'Created checkpoint directory: {ckpt_dir}')

        loader = self.load_dataset(otu_table, batch_size=batch_size)
        model = VAE(input_dim=otu_table.shape[1], latent_dim=model_configs['latent_dim']).to(device)
        # model = torch.nn.parallel.DataParallel(model)  # 使用DataParallel进行多GPU训练
        optimizer = torch.optim.Adam(model.parameters(), lr=lr)
        scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=epochs)  
        logger.info('Using CosineAnnealingLR scheduler with T_max set to the number of epochs.')

        logger.info('Starting VAE training...')
        epoch_losses: list = []
        epoch_recon_losses: list = []
        epoch_kl_losses: list = []
        epoch_lrs: list = []
        for epoch in range(1, epochs + 1):
            model.train()
            optimizer.zero_grad()
            step_losses: list = []
            step_recon_losses: list = []
            step_kl_losses: list = []
            for X_batch in loader:
                X_batch = X_batch[0].to(device)  # DataLoader returns a tuple, we need the first element
                recon_x, mu, logvar = model(X_batch)
                loss, recon_loss, kl = self.loss_function(recon_x, X_batch, mu, logvar, beta=model_configs['beta'])
                loss.backward()
                optimizer.step()
                optimizer.zero_grad()

                step_losses.append(loss.item())
                step_recon_losses.append(recon_loss.item())
                step_kl_losses.append(kl.item())

            epoch_loss = np.mean(step_losses)
            epoch_recon_loss = np.mean(step_recon_losses)
            epoch_kl_loss = np.mean(step_kl_losses)

            epoch_losses.append(epoch_loss)
            epoch_recon_losses.append(epoch_recon_loss)
            epoch_kl_losses.append(epoch_kl_loss)
            epoch_lrs.append(optimizer.param_groups[0]['lr'])

            scheduler.step()

            if epoch % 10 == 0 or epoch == 1 or epoch == epochs:
                print(
                    f'Epoch {epoch}, '
                    f'Train loss: {epoch_loss:.4f}, '
                    f'Recon loss: {epoch_recon_loss:.4f}, '
                    f'KL loss: {epoch_kl_loss:.4f}, '
                    f'lr: {optimizer.param_groups[0]["lr"]:.4e}')
                ckpt_dir_epoch = ckpt_dir.joinpath(f'vae_epoch_{epoch}.pth')
                torch.save(model.state_dict(), ckpt_dir_epoch)
            
        metric_df = pd.DataFrame({
            'Epoch': range(1, epochs + 1),
            'Train Loss': epoch_losses,
            'Recon Loss': epoch_recon_losses,
            'KL Loss': epoch_kl_losses,
            'Learning Rate': epoch_lrs,
        })
        logger.info('VAE training completed.')

        return model, metric_df
    
    def get_latent_representation(self, model: VAE, otu_table: pd.DataFrame, device: str = 'cpu') -> pd.DataFrame:
        '''Get the latent representation of the OTU table using the trained VAE model.'''

        model = model.to(device)
        model.eval()
        with torch.no_grad():
            X = torch.tensor(otu_table.values, dtype=torch.float32).to(device)
            recon_x, mu, logvar = model(X)
        
        latent_df = pd.DataFrame(data=mu.cpu().numpy(), index=otu_table.index)

        import umap
        umap_mu = umap.UMAP(n_components=2, random_state=42).fit_transform(mu.cpu().numpy())
        latent_df['UMAP1'] = umap_mu[:, 0]
        latent_df['UMAP2'] = umap_mu[:, 1]

        from sklearn.manifold import TSNE
        tsne_mu = TSNE(n_components=2, random_state=42).fit_transform(mu.cpu().numpy())
        latent_df['TSNE1'] = tsne_mu[:, 0]
        latent_df['TSNE2'] = tsne_mu[:, 1]
        
        from sklearn.decomposition import PCA
        pca_mu = PCA(n_components=2, random_state=42).fit_transform(mu.cpu().numpy())
        latent_df['PCA1'] = pca_mu[:, 0]
        latent_df['PCA2'] = pca_mu[:, 1]

        return latent_df
    
    def ecotype_clustering(self, latent_df: pd.DataFrame, method: str = 'kmeans', n_clusters: int = 4) -> pd.DataFrame:
        '''Cluster the latent representations into ecotypes using KMeans.'''

        new_latent_df = latent_df.copy()
        if method == 'kmeans':
            from sklearn.cluster import KMeans
            kmeans = KMeans(n_clusters=n_clusters, random_state=42)
            new_latent_df['Ecotype'] = kmeans.fit_predict(new_latent_df.drop(columns=['UMAP1', 'UMAP2', 'TSNE1', 'TSNE2', 'PCA1', 'PCA2']))
        elif method == 'gauss':
            from sklearn.mixture import GaussianMixture
            gmm = GaussianMixture(n_components=n_clusters, random_state=42)
            new_latent_df['Ecotype'] = gmm.fit_predict(new_latent_df.drop(columns=['UMAP1', 'UMAP2', 'TSNE1', 'TSNE2', 'PCA1', 'PCA2']))
        else:
            raise ValueError(f'Unsupported clustering method: {method}')

        return new_latent_df


class DownstreamTasks:
    '''Downstream tasks for evaluating the latent representations from the VAE.'''

    @staticmethod
    def classification(latent_df: pd.DataFrame, labels: pd.Series, models: list = ['rf', 'lr'], verbose: bool = True):
        '''Perform classification using the latent representations.'''

        from sklearn.model_selection import train_test_split, KFold, StratifiedKFold

        from sklearn.ensemble import RandomForestClassifier
        from sklearn.linear_model import LogisticRegression
        from sklearn.svm import SVC

        from sklearn.metrics import accuracy_score, precision_score, recall_score, f1_score
        from sklearn.metrics import classification_report, confusion_matrix


        kfold = StratifiedKFold(n_splits=10, shuffle=True, random_state=42).split(latent_df, labels)
        metrics: dict = {
            'model': [],
            'kfold': [],
            'accuracy': [],
            'precision': [],
            'recall': [],
            'f1': [],
        }
        for k, (train_idx, test_idx) in enumerate(kfold):
            X_train, X_test = latent_df.iloc[train_idx], latent_df.iloc[test_idx]
            y_train, y_test = labels.iloc[train_idx], labels.iloc[test_idx]

            for model in models:
                if model == 'rf':
                    clf = RandomForestClassifier(random_state=42)
                elif model == 'lr':
                    clf = LogisticRegression(random_state=42)
                elif model == 'svm':
                    clf = SVC(random_state=42)
                else:
                    raise ValueError(f'Unsupported model: {model}')
                clf.fit(X_train, y_train)
                y_pred = clf.predict(X_test)

                # report = classification_report(y_test, y_pred)
                # cm = confusion_matrix(y_test, y_pred)
                metrics['model'].append(model)
                metrics['kfold'].append(k+1)
                metrics['accuracy'].append(accuracy_score(y_test, y_pred))
                metrics['precision'].append(precision_score(y_test, y_pred, average='weighted'))
                metrics['recall'].append(recall_score(y_test, y_pred, average='weighted'))
                metrics['f1'].append(f1_score(y_test, y_pred, average='weighted'))

                if verbose:
                    print(
                        f'K-Fold {k+1}, Model: {model}:\n'
                        f'Accuracy: {metrics["accuracy"][-1]:.4f}, '
                        f'Precision: {metrics["precision"][-1]:.4f}, '
                        f'Recall: {metrics["recall"][-1]:.4f}, '
                        f'F1 Score: {metrics["f1"][-1]:.4f}'
                    )

        metrics_df = pd.DataFrame(metrics)
        print('\nOverall Classification Metrics:')
        print(metrics_df.describe().loc[['mean', 'std']])

        return metrics_df 

    @staticmethod
    def regression(latent_df: pd.DataFrame, targets: pd.Series, models: list = ['rf', 'lr'], verbose: bool = True):
        '''Perform regression using the latent representations.'''
        from sklearn.model_selection import train_test_split, KFold, StratifiedKFold
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.linear_model import LinearRegression
        from sklearn.svm import SVR
        from sklearn.metrics import mean_squared_error, r2_score

        kfold = KFold(n_splits=10, shuffle=True, random_state=42).split(latent_df)
        metrics: dict = {
            'model': [],
            'kfold': [],
            'mse': [],
            'r2': [],
        }
        for k, (train_idx, test_idx) in enumerate(kfold):
            X_train, X_test = latent_df.iloc[train_idx], latent_df.iloc[test_idx]
            y_train, y_test = targets.iloc[train_idx], targets.iloc[test_idx]

            for model in models:
                if model == 'rf':
                    reg = RandomForestRegressor(random_state=42)
                elif model == 'lr':
                    reg = LinearRegression()
                elif model == 'svm':
                    reg = SVR()
                else:
                    raise ValueError(f'Unsupported model: {model}')
                reg.fit(X_train, y_train)
                y_pred = reg.predict(X_test)

                metrics['model'].append(model)
                metrics['kfold'].append(k+1)
                metrics['mse'].append(mean_squared_error(y_test, y_pred))
                metrics['r2'].append(r2_score(y_test, y_pred))

                if verbose:
                    print(
                        f'K-Fold {k+1}, Model: {model}:\n'
                        f'MSE: {metrics["mse"][-1]:.4f}, '
                        f'R2 Score: {metrics["r2"][-1]:.4f}'
                    )
        metrics_df = pd.DataFrame(metrics)
        print('\nOverall Regression Metrics:')
        print(metrics_df.describe().loc[['mean', 'std']])

        return metrics_df

    @staticmethod
    def clustering(latent_df: pd.DataFrame, labels: pd.Series, method: str = 'kmeans', n_clusters: int = 4):
        '''Perform clustering using the latent representations and evaluate against true labels.'''

        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

        if method == 'kmeans':
            from sklearn.cluster import KMeans
            clusterer = KMeans(n_clusters=n_clusters, random_state=42)
            cluster_labels = clusterer.fit_predict(latent_df)
        elif method == 'gauss':
            from sklearn.mixture import GaussianMixture
            gmm = GaussianMixture(n_components=n_clusters, random_state=42)
            cluster_labels = gmm.fit_predict(latent_df)
        else:
            raise ValueError(f'Unsupported clustering method: {method}')

        ari = adjusted_rand_score(labels, cluster_labels)
        nmi = normalized_mutual_info_score(labels, cluster_labels)

        print(f'Clustering Evaluation - ARI: {ari:.4f}, NMI: {nmi:.4f}')

        return cluster_labels

    @staticmethod
    def beta_diversity(latent_df: pd.DataFrame, labels: pd.Series, metric: str = 'euclidean'):
        '''Calculate beta diversity metrics using the latent representations.'''

        from sklearn.metrics import pairwise_distances

        distance_matrix = pairwise_distances(latent_df, metric=metric)
        # 这里可以根据标签计算组内和组间的距离分布，或者使用 PERMANOVA 等方法进行统计检验

        return distance_matrix
    