"""
06_ols_mgwr.py
Actividad 15: regresión lineal múltiple (OLS) y, si los residuos presentan
autocorrelación espacial significativa, regresión geográficamente ponderada.

Para cada período (Ver24, Inv24, Ver25, Inv25), sobre la malla de 100 m:

  1. OLS  LST ~ NDVI + NDBI   (modelo de la metodología), más los modelos simples
     LST ~ NDVI y LST ~ NDBI para comparar por AICc (NDVI y NDBI son muy colineales).
     Reporta: coeficientes, p, R², R² ajustado, AICc, VIF y Moran's I de los
     residuos (Queen, estandarizada por filas, 999 permutaciones).

  2. Modelo local, solo si el Moran's I de los residuos del OLS es significativo:
       MODELO_LOCAL = "MGWR"  (recomendado) GWR multiescala: cada variable tiene su
                              propio ancho de banda (Fotheringham et al., 2017).
       MODELO_LOCAL = "GWR"   GWR clásico con un solo ancho de banda.
     Kernel bisquare adaptativo, ancho(s) de banda por AICc, variables estandarizadas.
     Reporta: anchos de banda, R², R² ajustado, AICc vs. OLS, Moran's I de los
     residuos, % de celdas con coeficiente significativo (t corregido por pruebas
     múltiples, da Silva y Fotheringham, 2016) y número de condición local.

Por qué MGWR: con el GWR clásico el ancho de banda óptimo cae al mínimo posible
(~15–47 vecinos) y los coeficientes cambian de signo de una celda a otra, porque un
solo ancho tiene que capturar a la vez la variación muy local de la LST y las
relaciones con NDVI/NDBI. En el MGWR el intercepto absorbe la variación local y
los coeficientes de NDVI y NDBI quedan estables e interpretables.

Tiempo aproximado: MGWR ~5–10 min por período (mgwr usa todos los núcleos del equipo).

Entradas: grid_lst_periodos.shp y ndvi_ndbi_malla100m.csv (los mismos de 05_correlaciones.py).
Salidas (carpeta SALIDA):
  ols_diagnosticos.csv, ols_coeficientes.csv, local_diagnosticos.csv,
  malla_modelo_local.gpkg (coeficientes, t, significancia y residuos por celda),
  mapas_<modelo>_<período>.png
"""

import warnings
from pathlib import Path

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from esda.moran import Moran
from libpysal.weights import Queen
from mgwr.gwr import GWR, MGWR
from mgwr.sel_bw import Sel_BW
from statsmodels.stats.outliers_influence import variance_inflation_factor

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================ PARÁMETROS ============================
# Carpeta raíz de los datos: por defecto, la carpeta que contiene a la de este script
BASE = Path(__file__).resolve().parents[1]
MALLA = BASE / "Clip_Grid" / "grid_lst_periodos.shp"
CSV_INDICES = BASE / "GEE_ICU" / "ndvi_ndbi_malla100m.csv"
SALIDA = BASE / "OLS_GWR"

MODELO_LOCAL = "MGWR"   # "MGWR" o "GWR"
CAMPO_ID = "id"
ALFA = 0.05
PERMUTACIONES = 999
SEMILLA = 12345

PERIODOS = {  # nombre corto: columna de LST en la malla
    "Ver24": "Ver24_mean",
    "Inv24": "Inv24_mean",
    "Ver25": "Ver25_mean",
    "Inv25": "Inv25_mean",
}
NOMBRES_LARGOS = {
    "Ver24": "Verano 2024", "Inv24": "Invierno 2024",
    "Ver25": "Verano 2025", "Inv25": "Invierno 2025",
}
VARIABLES = ["NDVI", "NDBI"]


# ============================ FUNCIONES ============================
def aicc_ols(modelo):
    """AICc del OLS (k = coeficientes + varianza del error), mismo criterio que mgwr."""
    n = int(modelo.nobs)
    k = int(modelo.df_model) + 2
    return modelo.aic + (2 * k * (k + 1)) / (n - k - 1)


def moran_residuos(residuos, pesos):
    np.random.seed(SEMILLA)
    mi = Moran(residuos, pesos, permutations=PERMUTACIONES)
    return mi.I, mi.z_sim, mi.p_sim


def ajustar_ols(df, variables, pesos):
    X = sm.add_constant(df[variables])
    modelo = sm.OLS(df["LST"], X).fit()
    I, z, p = moran_residuos(modelo.resid.values, pesos)
    if len(variables) > 1:
        vif = {v: variance_inflation_factor(X.values, i) for i, v in enumerate(variables, start=1)}
    else:
        vif = {variables[0]: 1.0}
    return modelo, {"R2": modelo.rsquared, "R2_ajustado": modelo.rsquared_adj,
                    "AICc": aicc_ols(modelo), "Moran_I": I, "Moran_z": z, "Moran_p": p, "vif": vif}


def ajustar_local(coords, ys, Xs):
    """Ajusta MGWR o GWR con variables estandarizadas. Devuelve resultados y anchos de banda."""
    if MODELO_LOCAL == "MGWR":
        selector = Sel_BW(coords, ys, Xs, multi=True, kernel="bisquare", fixed=False)
        anchos = selector.search(multi_bw_min=[2])
        res = MGWR(coords, ys, Xs, selector, kernel="bisquare", fixed=False).fit()
        cn = res.local_collinearity()[0].flatten()
        r2_local = None
    else:
        selector = Sel_BW(coords, ys, Xs, kernel="bisquare", fixed=False)
        ancho = selector.search(criterion="AICc")
        anchos = np.repeat(ancho, Xs.shape[1] + 1)
        res = GWR(coords, ys, Xs, bw=ancho, kernel="bisquare", fixed=False).fit()
        cn = res.local_collinearity()[2].flatten()
        r2_local = res.localR2.flatten()
    return res, np.asarray(anchos).flatten(), cn, r2_local


def mapa_local(malla, periodo, archivo):
    """Coeficientes locales estandarizados (gris = no significativos)."""
    columnas = ["b0"] + VARIABLES
    titulos = {"b0": "Intercepto local", "NDVI": "Coeficiente local NDVI", "NDBI": "Coeficiente local NDBI"}
    fig, ejes = plt.subplots(1, 3, figsize=(15, 7.5), constrained_layout=True)
    for eje, var in zip(ejes, columnas):
        col = f"{var}_std_{periodo}"
        datos = malla[[col, f"sig_{var}_{periodo}", "geometry"]].copy()
        datos.loc[datos[f"sig_{var}_{periodo}"] == 0, col] = np.nan
        lim = np.nanmax(np.abs(datos[col])) if datos[col].notna().any() else 1
        datos.plot(column=col, cmap="RdBu_r", vmin=-lim, vmax=lim, ax=eje, legend=True,
                   legend_kwds={"shrink": 0.6, "label": "Coeficiente estandarizado"},
                   missing_kwds={"color": "#d9d9d9"})
        eje.set_title(titulos[var], fontsize=11)
        eje.set_axis_off()
    fig.suptitle(f"{MODELO_LOCAL} – LST ~ NDVI + NDBI – {NOMBRES_LARGOS[periodo]}\n"
                 f"Gris: coeficiente no significativo (t corregido, α = {ALFA})", fontsize=12)
    fig.savefig(archivo, dpi=300)
    plt.close(fig)


# ============================ EJECUCIÓN ============================
def main():
    SALIDA.mkdir(parents=True, exist_ok=True)

    malla = gpd.read_file(MALLA)
    malla = malla.merge(pd.read_csv(CSV_INDICES), on=CAMPO_ID, how="inner").reset_index(drop=True)
    print(f"Celdas en el análisis: {len(malla)}  |  Modelo local: {MODELO_LOCAL}")

    pesos = Queen.from_dataframe(malla, use_index=True)
    pesos.transform = "r"
    coords = np.column_stack([malla.geometry.centroid.x, malla.geometry.centroid.y])

    diag_ols, coef_ols, diag_local = [], [], []

    for periodo, col_lst in PERIODOS.items():
        print(f"\n===== {NOMBRES_LARGOS[periodo]} =====")
        df = pd.DataFrame({"LST": malla[col_lst],
                           "NDVI": malla[f"NDVI_{periodo}"],
                           "NDBI": malla[f"NDBI_{periodo}"]})
        if df.isna().any().any():
            raise SystemExit(f"Hay valores vacíos en {periodo}; revisá la malla o el CSV.")

        # ---------------- OLS ----------------
        modelos = {"NDVI + NDBI": VARIABLES, "Solo NDVI": ["NDVI"], "Solo NDBI": ["NDBI"]}
        for nombre, variables in modelos.items():
            modelo, d = ajustar_ols(df, variables, pesos)
            if nombre == "NDVI + NDBI":
                ols, d_ols = modelo, d
            diag_ols.append({
                "periodo": periodo, "modelo": nombre, "n": int(modelo.nobs),
                "R2": d["R2"], "R2_ajustado": d["R2_ajustado"], "AICc": d["AICc"],
                "VIF_NDVI": d["vif"].get("NDVI", np.nan), "VIF_NDBI": d["vif"].get("NDBI", np.nan),
                "Moran_I_residuos": d["Moran_I"], "Moran_z": d["Moran_z"], "Moran_p": d["Moran_p"],
            })
            for var in ["const"] + variables:
                coef_ols.append({"periodo": periodo, "modelo": nombre, "variable": var,
                                 "coeficiente": modelo.params[var], "error_estandar": modelo.bse[var],
                                 "t": modelo.tvalues[var], "p": modelo.pvalues[var]})
            print(f"OLS {nombre:12s} R²={d['R2']:.3f}  R²aj={d['R2_ajustado']:.3f}  "
                  f"AICc={d['AICc']:.1f}  Moran I res={d['Moran_I']:.3f} (p={d['Moran_p']:.3f})")
        print(f"VIF (NDVI, NDBI): {d_ols['vif']['NDVI']:.2f}, {d_ols['vif']['NDBI']:.2f}")

        # ---------------- Modelo local (condicionado) ----------------
        if d_ols["Moran_p"] >= ALFA:
            print("Residuos del OLS sin autocorrelación significativa: no se aplica modelo local.")
            continue

        y = df[["LST"]].values
        X = df[VARIABLES].values
        ys = (y - y.mean()) / y.std()
        Xs = (X - X.mean(axis=0)) / X.std(axis=0)

        # Para comparar AICc con el OLS en la misma escala, OLS estandarizado
        ols_std = sm.OLS(ys, sm.add_constant(Xs)).fit()
        aicc_ols_std = aicc_ols(ols_std)

        print(f"Ajustando {MODELO_LOCAL} (puede tardar varios minutos)...")
        res, anchos, cn, r2_local = ajustar_local(coords, ys, Xs)
        I_l, z_l, p_l = moran_residuos(res.resid_response.flatten(), pesos)
        filtrados = res.filter_tvals(alpha=ALFA)

        fila = {
            "periodo": periodo, "modelo": MODELO_LOCAL, "n": res.n,
            "ancho_intercepto": int(anchos[0]), "ancho_NDVI": int(anchos[1]), "ancho_NDBI": int(anchos[2]),
            "R2": res.R2, "R2_ajustado": res.adj_R2,
            "AICc": res.aicc, "AICc_OLS_estandarizado": aicc_ols_std,
            "dif_AICc_OLS_menos_local": aicc_ols_std - res.aicc,
            "Moran_I_residuos": I_l, "Moran_z": z_l, "Moran_p": p_l,
            "pct_celdas_numero_condicion_mayor_30": 100 * np.mean(cn > 30),
        }
        sd_y = y.std()
        for j, var in enumerate(["b0"] + VARIABLES):
            b = res.params[:, j]
            fila[f"{var}_std_mediana"] = np.median(b)
            fila[f"{var}_std_p05"] = np.percentile(b, 5)
            fila[f"{var}_std_p95"] = np.percentile(b, 95)
            fila[f"{var}_pct_positivo"] = 100 * np.mean(b > 0)
            fila[f"{var}_pct_signif"] = 100 * np.mean(filtrados[:, j] != 0)
            malla[f"{var}_std_{periodo}"] = b
            malla[f"t_{var}_{periodo}"] = res.tvalues[:, j]
            malla[f"sig_{var}_{periodo}"] = (filtrados[:, j] != 0).astype(int)
            if var != "b0":
                # Coeficiente en unidades originales: °C por unidad de índice
                malla[f"{var}_degC_{periodo}"] = b * sd_y / X[:, j - 1].std()
                fila[f"{var}_degC_mediana"] = np.median(malla[f"{var}_degC_{periodo}"])
        diag_local.append(fila)

        malla[f"res_local_{periodo}"] = res.resid_response.flatten() * sd_y  # en °C
        malla[f"res_OLS_{periodo}"] = ols.resid.values
        malla[f"CN_{periodo}"] = cn
        if r2_local is not None:
            malla[f"R2loc_{periodo}"] = r2_local

        print(f"{MODELO_LOCAL} anchos (intercepto, NDVI, NDBI) = {anchos.astype(int).tolist()}  "
              f"R²={res.R2:.3f}  R²aj={res.adj_R2:.3f}  AICc={res.aicc:.1f} "
              f"(OLS {aicc_ols_std:.1f})  Moran I res={I_l:.3f} (p={p_l:.3f})")

        mapa_local(malla, periodo, SALIDA / f"mapas_{MODELO_LOCAL}_{periodo}.png")

    pd.DataFrame(diag_ols).to_csv(SALIDA / "ols_diagnosticos.csv", index=False)
    pd.DataFrame(coef_ols).to_csv(SALIDA / "ols_coeficientes.csv", index=False)
    if diag_local:
        pd.DataFrame(diag_local).to_csv(SALIDA / "local_diagnosticos.csv", index=False)
    malla.to_file(SALIDA / "malla_modelo_local.gpkg", driver="GPKG", layer=f"{MODELO_LOCAL}_resultados")
    print(f"\nListo. Resultados en: {SALIDA}")


if __name__ == "__main__":
    main()
