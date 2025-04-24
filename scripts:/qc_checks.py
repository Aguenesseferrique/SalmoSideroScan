import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import os

# === Paramètres ===
# Remplacez ce chemin par le chemin local de votre fichier CSV
input_csv = '20240307_Salmonella_genomes.csv'
output_ids_file = 'GOODSALMONELLA1.txt'
plot_png = 'CDSLinear_Contam.png'
plot_svg = 'CDSLinear_Contam.svg'

# === Traitement ===
try:
    # Lecture du fichier CSV
    if not os.path.isfile(input_csv):
        raise FileNotFoundError(f"Fichier non trouvé : {input_csv}")
    
    df = pd.read_csv(input_csv, sep=';', header=0)
    print(f"Colonnes détectées : {df.columns.tolist()}")

    # Colonnes de métadonnées à vérifier
    metadata_columns = ['Isolation Source', 'Isolation Comments', 'Host Name', 'Host Common Name']
    missing_columns = [col for col in metadata_columns if col not in df.columns]

    if missing_columns:
        raise ValueError(f"Colonnes manquantes dans le fichier : {', '.join(missing_columns)}")

    # Filtrage des lignes avec métadonnées absentes
    df = df.dropna(subset=metadata_columns, how='all')

    # Filtrage CheckM
    df = df.loc[(df['CheckM Completeness'] > 90) | df['CheckM Completeness'].isna()]
    df = df.loc[(df['CheckM Contamination'] < 5) | df['CheckM Contamination'].isna()]

    # Régression : longueur vs CDS
    df['lengthkbp'] = df['Size'] / 1000
    m, b = np.polyfit(x=df['lengthkbp'], y=df['CDS'], deg=1)
    df['predicted'] = m * df['lengthkbp'] + b
    df['error'] = abs(df['CDS'] - df['predicted'])
    df['Remove'] = np.where(df['error'] > 500, 'Yes', 'No')

    # Plot
    sns.set(style='white', font_scale=1.6, palette='pastel')
    sns.lmplot(x='lengthkbp', y='CDS', hue='Remove', fit_reg=False, data=df,
               scatter_kws={'edgecolors': 'k', 's': 80})
    ax = sns.regplot(x='lengthkbp', y='CDS', data=df, scatter_kws={"zorder": -1},
                     line_kws={'label': f"y={m:.3f}x+{b:.1f}", 'color': 'k'})
    ax.figure.set_size_inches(18.5, 10.5)
    ax.legend()
    plt.ylabel('Coding Sequences')
    plt.xlabel('Genome Length (kbp)')
    plt.savefig(plot_png)
    plt.savefig(plot_svg)

    # Export des IDs valides
    df.loc[df['Remove'] == 'No', 'Genome ID'].to_csv(output_ids_file, index=False, header=False)
    print(f"IDs valides sauvegardés dans {output_ids_file}")

except pd.errors.ParserError as e:
    print(f"Erreur de parsing du CSV : {e}")
except ValueError as e:
    print(f"Erreur de validation : {e}")
except FileNotFoundError as e:
    print(f"Erreur de fichier : {e}")
except Exception as e:
    print(f"Erreur inattendue : {e}")
