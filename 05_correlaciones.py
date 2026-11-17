"""
05_correlaciones.py
Actividad 14: correlaciones de Pearson y Spearman entre la LST y las variables
ambientales, en las dos escalas definidas en la metodología (5.1.3):

  A) Nivel de celda (malla de 100 m), una matriz por período:
     LST, NDVI y NDBI.
  B) Nivel zona × mes, una matriz por período:
     LST, NDVI, NDBI, NO2, CO y PM2.5 (cumple el mínimo de 4 variables
     del indicador del Objetivo 2).
     Además: correlación temporal LST–contaminantes por zona (2024–2025).
  C) Tabla de contexto meteorológico (ERA5-Land) por período.

Entradas:
  - grid_lst_periodos.shp  (LST por período en cada celda)
  - ndvi_ndbi_malla100m.csv  (salida de gee/04_ndvi_ndbi_malla100m.js)
  - serie_mensual_zonas_2024_2025.csv  (salida de gee/03_serie_mensual_zonas.js)

Salidas (carpeta SALIDA): tablas CSV con r, p y n, y figuras PNG a 300 dpi.
"""

from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

# ============================ PARÁMETROS ============================
# Carpeta raíz de los datos: por defecto, la carpeta que contiene a la de este script
BASE = Path(__file__).resolve().parents[1]
MALLA = BASE / "Clip_Grid" / "grid_lst_periodos.shp"
CSV_INDICES = BASE / "GEE_ICU" / "ndvi_ndbi_malla100m.csv"
CSV_SERIE = BASE / "GEE_ICU" / "serie_mensual_zonas_2024_2025.csv"
SALIDA = BASE / "Correlaciones"

CAMPO_ID = "id"  # identificador de celda, igual que en el script 02
ALFA = 0.05
UMBRAL_LST_VALIDO = 0.5  # fracción mínima de pixeles sin nubes para aceptar la LST mensual de una zona

# Nombre de la columna de LST de cada período dentro de grid_lst_periodos.shp.
# El script imprime las columnas disponibles: ajustá estos nombres si no coinciden.
COLUMNAS_LST = {
    "Ver24": "Ver24_mean",
    "Inv24": "Inv24_mean",
    "Ver25": "Ver25_mean",
    "Inv25": "Inv25_mean",
}

# Correspondencia entre los nombres cortos y los del CSV mensual
PERIODOS_SERIE = {
    "Ver24": "Verano2024",
    "Inv24": "Invierno2024",
    "Ver25": "Verano2025",
    "Inv25": "Invierno2025",
}

ETIQUETAS = {
    "LST": "LST (°C)",
    "NDVI": "NDVI",
    "NDBI": "NDBI",
    "NO2": "NO₂ (µmol/m²)",
    "CO": "CO (mmol/m²)",
    "PM25": "PM₂.₅ (µg/m³)",
}


# ============================ FUNCIONES ============================
def matriz_correlacion(df, variables, metodo):
    """Devuelve tablas de r, p y n para todas las parejas de variables."""
    r = pd.DataFrame(np.nan, index=variables, columns=variables)
    p = r.copy()
    n = pd.DataFrame(0, index=variables, columns=variables)
    for i, a in enumerate(variables):
        for b in variables[i:]:
            par = df[[a, b]].dropna()
            n.loc[a, b] = n.loc[b, a] = len(par)
            if a == b:
                r.loc[a, b], p.loc[a, b] = 1.0, 0.0
                continue
            if len(par) < 3 or par[a].nunique() < 2 or par[b].nunique() < 2:
                continue
            if metodo == "pearson":
                res = stats.pearsonr(par[a], par[b])
            else:
                res = stats.spearmanr(par[a], par[b])
            r.loc[a, b] = r.loc[b, a] = res[0]
            p.loc[a, b] = p.loc[b, a] = res[1]
    return r, p, n


def tabla_larga(r, p, n, periodo, metodo, escala):
    """Convierte las matrices en una tabla larga (una fila por pareja)."""
    filas = []
    vars_ = list(r.index)
    for i, a in enumerate(vars_):
        for b in vars_[i + 1:]:
            filas.append({
                "escala": escala,
                "periodo": periodo,
                "metodo": metodo,
                "variable_1": a,
                "variable_2": b,
                "r": r.loc[a, b],
                "p": p.loc[a, b],
                "n": n.loc[a, b],
                "significativo": bool(p.loc[a, b] < ALFA) if pd.notna(p.loc[a, b]) else False,
            })
    return filas


def figura_matrices(resultados, titulo, archivo):
    """Panel 2×2 de matrices (una por período) con r y asterisco si p < α."""
    fig, ejes = plt.subplots(2, 2, figsize=(11, 10), constrained_layout=True)
    for eje, (periodo, (r, p)) in zip(ejes.flat, resultados.items()):
        vals = r.values.astype(float)
        im = eje.imshow(vals, cmap="RdBu_r", vmin=-1, vmax=1)
        etiquetas = [ETIQUETAS.get(v, v) for v in r.index]
        eje.set_xticks(range(len(etiquetas)), etiquetas, rotation=45, ha="right")
        eje.set_yticks(range(len(etiquetas)), etiquetas)
        for i in range(len(etiquetas)):
            for j in range(len(etiquetas)):
                if np.isnan(vals[i, j]):
                    texto = "–"
                else:
                    sig = "*" if (i != j and p.values[i, j] < ALFA) else ""
                    texto = f"{vals[i, j]:.2f}{sig}"
                color = "white" if (not np.isnan(vals[i, j]) and abs(vals[i, j]) > 0.6) else "black"
                eje.text(j, i, texto, ha="center", va="center", fontsize=9, color=color)
        eje.set_title(periodo)
    fig.colorbar(im, ax=ejes, shrink=0.6, label="Coeficiente de correlación")
    fig.suptitle(f"{titulo}\n* p < {ALFA}", fontsize=12)
    fig.savefig(archivo, dpi=300)
    plt.close(fig)


# ============================ A) NIVEL DE CELDA ============================
def correlacion_celdas():
    malla = gpd.read_file(MALLA)
    print("Columnas en la malla:", list(malla.columns))
    faltan = [c for c in list(COLUMNAS_LST.values()) + [CAMPO_ID] if c not in malla.columns]
    if faltan:
        raise SystemExit(f"Ajustá COLUMNAS_LST / CAMPO_ID; no encuentro: {faltan}")

    indices = pd.read_csv(CSV_INDICES)
    datos = malla.drop(columns="geometry").merge(indices, on=CAMPO_ID, how="inner")
    print(f"Celdas unidas: {len(datos)} de {len(malla)}")

    filas, figs = [], {"pearson": {}, "spearman": {}}
    for corto, col_lst in COLUMNAS_LST.items():
        df = pd.DataFrame({
            "LST": datos[col_lst],
            "NDVI": datos[f"NDVI_{corto}"],
            "NDBI": datos[f"NDBI_{corto}"],
        })
        for metodo in ("pearson", "spearman"):
            r, p, n = matriz_correlacion(df, ["LST", "NDVI", "NDBI"], metodo)
            filas += tabla_larga(r, p, n, corto, metodo, "celda_100m")
            figs[metodo][corto] = (r, p)

        # Dispersión LST vs NDVI y NDBI
        fig, ejes = plt.subplots(1, 2, figsize=(10, 4.2), constrained_layout=True)
        for eje, var in zip(ejes, ("NDVI", "NDBI")):
            par = df[["LST", var]].dropna()
            eje.scatter(par[var], par["LST"], s=4, alpha=0.35, color="#3A4750")
            if len(par) > 2:
                pend, inter, rr, pp, _ = stats.linregress(par[var], par["LST"])
                xs = np.linspace(par[var].min(), par[var].max(), 50)
                eje.plot(xs, inter + pend * xs, color="#d73027", lw=1.5)
                eje.set_title(f"r = {rr:.2f}, n = {len(par)}")
            eje.set_xlabel(var)
            eje.set_ylabel("LST (°C)")
        fig.suptitle(f"LST vs. índices espectrales – {corto}")
        fig.savefig(SALIDA / f"dispersion_celda_{corto}.png", dpi=300)
        plt.close(fig)

    for metodo in ("pearson", "spearman"):
        figura_matrices(figs[metodo], f"Correlación de {metodo.capitalize()} – malla de 100 m",
                        SALIDA / f"matrices_celda_{metodo}.png")
    return filas


# ============================ B) ZONA × MES ============================
def depurar_serie(serie):
    """Control de calidad: descarta LST con poca cobertura sin nubes y NDVI/NDBI no plausibles."""
    serie = serie.copy()
    serie["zona"] = serie["zona"].astype(str)
    if "pct_LST_valido" in serie.columns:
        malos = serie["pct_LST_valido"].fillna(0) < UMBRAL_LST_VALIDO
        serie.loc[malos, "LST"] = np.nan
        print(f"LST descartada por cobertura < {UMBRAL_LST_VALIDO:.0%}: {int(malos.sum())} filas zona-mes")
    # Un NDVI medio de zona urbana ≤ 0 indica contaminación por nubes (pocas escenas Sentinel-2)
    malos_ndvi = serie["NDVI"] <= 0
    serie.loc[malos_ndvi, ["NDVI", "NDBI"]] = np.nan
    print(f"NDVI/NDBI descartados por valor no plausible: {int(malos_ndvi.sum())} filas zona-mes")
    serie.to_csv(SALIDA / "serie_mensual_depurada.csv", index=False)
    return serie


def correlacion_zona_mes():
    serie = depurar_serie(pd.read_csv(CSV_SERIE))
    zonas = serie[serie["zona"] != "Total"].copy()
    variables = [v for v in ["LST", "NDVI", "NDBI", "NO2", "CO", "PM25"] if v in zonas.columns]

    filas, figs = [], {"pearson": {}, "spearman": {}}
    for corto, largo in PERIODOS_SERIE.items():
        df = zonas[zonas["periodo"] == largo]
        for metodo in ("pearson", "spearman"):
            r, p, n = matriz_correlacion(df, variables, metodo)
            filas += tabla_larga(r, p, n, corto, metodo, "zona_mes")
            figs[metodo][corto] = (r, p)

    for metodo in ("pearson", "spearman"):
        figura_matrices(figs[metodo], f"Correlación de {metodo.capitalize()} – zona × mes",
                        SALIDA / f"matrices_zona_mes_{metodo}.png")

    # Correlación temporal por zona (todos los meses de las épocas analizadas)
    estacional = zonas[~zonas["periodo"].str.startswith("Transicion")]
    temporales = []
    for zona, df in estacional.groupby("zona"):
        for contaminante in ("NO2", "CO", "PM25"):
            par = df[["LST", contaminante]].dropna()
            for metodo, fn in (("pearson", stats.pearsonr), ("spearman", stats.spearmanr)):
                if len(par) >= 3:
                    res = fn(par["LST"], par[contaminante])
                    r, p = res[0], res[1]
                else:
                    r, p = np.nan, np.nan
                temporales.append({"zona": zona, "contaminante": contaminante, "metodo": metodo,
                                   "r": r, "p": p, "n_meses": len(par),
                                   "significativo": bool(p < ALFA) if pd.notna(p) else False})
    pd.DataFrame(temporales).to_csv(SALIDA / "correlacion_temporal_por_zona.csv", index=False)

    # Series temporales: LST vs. contaminantes (área total)
    total = serie[serie["zona"] == "Total"].sort_values("fecha")
    if "LST" in total.columns:
        fig, ejes = plt.subplots(3, 1, figsize=(10, 8), sharex=True, constrained_layout=True)
        fechas = pd.to_datetime(total["fecha"])
        for eje, cont in zip(ejes, ("NO2", "CO", "PM25")):
            eje.plot(fechas, total["LST"], marker="o", color="#d73027", label="LST")
            eje.set_ylabel("LST (°C)", color="#d73027")
            eje2 = eje.twinx()
            eje2.plot(fechas, total[cont], marker="s", color="#2166ac", label=cont)
            eje2.set_ylabel(ETIQUETAS[cont], color="#2166ac")
        fig.suptitle("Serie mensual 2024–2025: LST y contaminantes (zonas 1, 4, 9 y 10)")
        fig.savefig(SALIDA / "serie_LST_contaminantes.png", dpi=300)
        plt.close(fig)
    return filas


# ============================ C) CONTEXTO ERA5 ============================
def contexto_meteorologico():
    serie = pd.read_csv(CSV_SERIE)
    serie["zona"] = serie["zona"].astype(str)
    total = serie[(serie["zona"] == "Total") & (serie["periodo"].isin(PERIODOS_SERIE.values()))]
    tabla = (total.groupby("periodo")[["T_aire", "HR", "Viento"]]
             .agg(["mean", "min", "max"]).round(2)
             .reindex(PERIODOS_SERIE.values()))
    tabla.to_csv(SALIDA / "contexto_meteorologico_ERA5.csv")
    print("\nContexto meteorológico por período (ERA5-Land):\n", tabla)


# ============================ EJECUCIÓN ============================
if __name__ == "__main__":
    SALIDA.mkdir(parents=True, exist_ok=True)
    filas = correlacion_celdas() + correlacion_zona_mes()
    resumen = pd.DataFrame(filas)
    resumen.to_csv(SALIDA / "correlaciones_resumen.csv", index=False)
    contexto_meteorologico()
    print(f"\nListo. Resultados en: {SALIDA}")
