"""
07_mapas_tematicos.py
Mapas temáticos finales (Fase IV) con la plantilla del proyecto:
Carta horizontal (279.4 × 215.9 mm), banda de título, mapa principal a 1:50 000 con
cuadrícula en grados-minutos-segundos, norte, escala gráfica y numérica, mapas de
ubicación nacional y municipal, leyenda y cuadro técnico.

Mapas que genera (PDF y PNG a 300 dpi):
  - ICU por intensidad (clases μ ± σ), 4 períodos
  - LST estacional, 4 períodos (misma escala de color en los cuatro)
  - NDVI estacional, 4 períodos
  - Coeficiente local NDVI del MGWR (relación LST–NDVI), 4 períodos
"""

import sys
from pathlib import Path

import geopandas as gpd
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import BoundaryNorm, ListedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from pyproj import Transformer
from shapely.geometry import box

# ============================ PARÁMETROS ============================
BASE = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parents[1]
SALIDA = Path(sys.argv[2]) if len(sys.argv) > 2 else BASE / "Mapas_finales"
PAISES = Path(sys.argv[3]) if len(sys.argv) > 3 else BASE / "naturalearth_lowres" / "naturalearth_lowres.shp"
SOLO = sys.argv[4].split(",") if len(sys.argv) > 4 else None   # p.ej. "ICU_Ver24" para generar uno solo

AUTOR = "Jorge Cordón, Universidad Rafael Landívar, 2026"
METODO_LST = "Split-Window · 30 m"   # ajustar si el método final es otro
ESCALA = 50_000
CRS = "EPSG:32615"

PERIODOS = {
    "Ver24": ("Verano 2024", "Verano2024", "Verano_2024", "enero–abril 2024 (verano)"),
    "Inv24": ("Invierno 2024", "Invierno2024", "Invierno_2024", "mayo–octubre 2024 (invierno)"),
    "Ver25": ("Verano 2025", "Verano2025", "Verano_2025", "enero–abril 2025 (verano)"),
    "Inv25": ("Invierno 2025", "Invierno2025", "Invierno_2025", "mayo–octubre 2025 (invierno)"),
}

TITULO_FONDO = "#3A4750"
TINTA = "#1f1f1f"
TINTA_2 = "#555555"
FUENTE = "DejaVu Sans"
MM = 1 / 25.4  # pulgadas por mm
ANCHO, ALTO = 279.4, 215.9

plt.rcParams.update({"font.family": FUENTE, "font.size": 7})


# ============================ UTILIDADES ============================
def rect(x, y, w, h):
    """Rectángulo en mm (origen arriba-izquierda, como en QGIS) → fracciones de figura."""
    return [x / ANCHO, 1 - (y + h) / ALTO, w / ANCHO, h / ALTO]


def leer_raster(ruta, banda=1):
    with rasterio.open(ruta) as r:
        a = r.read(banda).astype(float)
        if r.nodata is not None:
            a[a == r.nodata] = np.nan
        ext = [r.bounds.left, r.bounds.right, r.bounds.bottom, r.bounds.top]
    return a, ext


def dms(valor, eje):
    v = abs(valor)
    g = int(v)
    m = int((v - g) * 60)
    s = round(((v - g) * 60 - m) * 60)
    if s == 60:
        s, m = 0, m + 1
    if m == 60:
        m, g = 0, g + 1
    hem = ("N" if valor >= 0 else "S") if eje == "lat" else ("E" if valor >= 0 else "O")
    return f"{g}°{m:02d}'{s:02d}\" {hem}"


def cuadricula_dms(ax, ext, paso=0.0125):
    """Líneas de latitud/longitud (WGS84) sobre el mapa en UTM, con etiquetas fuera del marco."""
    a_geo = Transformer.from_crs(CRS, "EPSG:4326", always_xy=True)
    a_utm = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    x0, x1, y0, y1 = ext
    lons, lats = a_geo.transform([x0, x1, x0, x1], [y0, y0, y1, y1])
    lon_min, lon_max, lat_min, lat_max = min(lons), max(lons), min(lats), max(lats)
    estilo = dict(color="#9a9a9a", lw=0.35, zorder=8)
    lat_s = np.linspace(lat_min - 0.01, lat_max + 0.01, 200)
    for lon in np.arange(np.ceil(lon_min / paso) * paso, lon_max, paso):
        xs, ys = a_utm.transform(np.full_like(lat_s, lon), lat_s)
        ax.plot(xs, ys, **estilo)
        xb = np.interp(y0, ys, xs)
        xt = np.interp(y1, ys, xs)
        if x0 < xb < x1:
            ax.text(xb, y0, dms(lon, "lon"), ha="center", va="top", fontsize=5.5, color=TINTA_2,
                    transform=ax.transData, clip_on=False).set_y(y0 - (y1 - y0) * 0.008)
        if x0 < xt < x1:
            ax.text(xt, y1 + (y1 - y0) * 0.008, dms(lon, "lon"), ha="center", va="bottom", fontsize=5.5,
                    color=TINTA_2, clip_on=False)
    lon_s = np.linspace(lon_min - 0.01, lon_max + 0.01, 200)
    for lat in np.arange(np.ceil(lat_min / paso) * paso, lat_max, paso):
        xs, ys = a_utm.transform(lon_s, np.full_like(lon_s, lat))
        ax.plot(xs, ys, **estilo)
        yl = np.interp(x0, xs, ys)
        yr = np.interp(x1, xs, ys)
        if y0 < yl < y1:
            ax.text(x0 - (x1 - x0) * 0.01, yl, dms(lat, "lat"), ha="right", va="center", fontsize=5.5,
                    color=TINTA_2, rotation=90, clip_on=False)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)


def norte(ax, x=0.93, y=0.95):
    ax.annotate("", xy=(x, y), xytext=(x, y - 0.07), xycoords="axes fraction",
                arrowprops=dict(facecolor=TINTA, edgecolor=TINTA, width=3, headwidth=9, headlength=9),
                zorder=20)
    ax.text(x, y + 0.008, "N", transform=ax.transAxes, ha="center", va="bottom", fontsize=10,
            fontweight="bold", color=TINTA, zorder=20)


def escala_grafica(ax, ext, km=3):
    x0, x1, y0, y1 = ext
    bx = x0 + (x1 - x0) * 0.05
    by = y0 + (y1 - y0) * 0.045
    alto = (y1 - y0) * 0.008
    for i in range(km):
        ax.add_patch(mpatches.Rectangle((bx + i * 1000, by), 1000, alto,
                                        facecolor=TINTA if i % 2 == 0 else "white",
                                        edgecolor=TINTA, lw=0.6, zorder=20))
        ax.text(bx + i * 1000, by + alto * 1.6, f"{i}", ha="center", va="bottom", fontsize=6, zorder=20)
    ax.text(bx + km * 1000, by + alto * 1.6, f"{km} km", ha="center", va="bottom", fontsize=6, zorder=20)
    ax.text(bx, by - alto * 1.4, f"Escala 1:{ESCALA:,}".replace(",", " "), ha="left", va="top",
            fontsize=6.5, zorder=20)


def encuadre(ax, bounds, ancho_mm, alto_mm, margen=0.08):
    """Ajusta los límites para que el área quepa completa respetando la proporción del recuadro."""
    x0, y0, x1, y1 = bounds
    w, h = (x1 - x0) * (1 + margen), (y1 - y0) * (1 + margen)
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    r = ancho_mm / alto_mm
    if w / h > r:
        h = w / r
    else:
        w = h * r
    ax.set_xlim(cx - w / 2, cx + w / 2)
    ax.set_ylim(cy - h / 2, cy + h / 2)
    ax.set_aspect("auto")


def marco(ax, titulo=None):
    for s in ax.spines.values():
        s.set_edgecolor(TINTA)
        s.set_linewidth(0.6)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    if titulo:
        ax.set_title(titulo, fontsize=7, pad=2, color=TINTA, loc="left")


# ============================ CAPAS BASE ============================
zonas_est = gpd.read_file(BASE / "zonas14910" / "Zonas14910_Guate_city.shp").to_crs(CRS)
zonas_mun = gpd.read_file(BASE / "zonas" / "zonas.shp")
if zonas_mun.crs is None:
    zonas_mun = zonas_mun.set_crs(CRS)
zonas_mun = zonas_mun.to_crs(CRS)
paises = gpd.read_file(PAISES)

cx = (zonas_est.total_bounds[0] + zonas_est.total_bounds[2]) / 2
cy = (zonas_est.total_bounds[1] + zonas_est.total_bounds[3]) / 2
MAPA_W, MAPA_H = 159.0, 176.0
EXT = [cx - MAPA_W / 1000 * ESCALA / 2, cx + MAPA_W / 1000 * ESCALA / 2,
       cy - MAPA_H / 1000 * ESCALA / 2, cy + MAPA_H / 1000 * ESCALA / 2]


def lamina(titulo, dibujar, leyenda, periodo_txt, tecnico_extra, archivo):
    fig = plt.figure(figsize=(ANCHO * MM, ALTO * MM))
    fig.patch.set_facecolor("white")

    # Banda de título
    axt = fig.add_axes(rect(10, 8, 259.4, 14))
    axt.set_facecolor(TITULO_FONDO)
    axt.text(0.5, 0.5, titulo, ha="center", va="center", color="white", fontsize=13, fontweight="bold")
    marco(axt)
    for s in axt.spines.values():
        s.set_visible(False)

    # Mapa principal
    ax = fig.add_axes(rect(16, 30, MAPA_W, MAPA_H))
    ax.set_facecolor("white")
    zonas_mun.plot(ax=ax, facecolor="#f2f2f2", edgecolor="white", lw=0.8, zorder=1)
    dibujar(ax)
    zonas_est.plot(ax=ax, facecolor="none", edgecolor=TINTA, lw=1.0, zorder=10)
    for _, z in zonas_est.iterrows():
        p = z.geometry.representative_point()
        ax.text(p.x, p.y, f"Zona {z['ZONA']}", ha="center", va="center", fontsize=6.5, color=TINTA,
                fontweight="bold", zorder=11,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="none", alpha=0.75))
    cuadricula_dms(ax, EXT)
    ax.set_aspect("equal")
    ax.set_xlim(EXT[0], EXT[1])
    ax.set_ylim(EXT[2], EXT[3])
    marco(ax)
    norte(ax)
    escala_grafica(ax, EXT)

    # Ubicación nacional
    axn = fig.add_axes(rect(180, 26, 89.4, 50))
    gt = paises[paises["name"] == "Guatemala"]
    vecinos = paises[paises["name"].isin(["Mexico", "Belize", "Honduras", "El Salvador"])]
    vecinos.plot(ax=axn, facecolor="#ededed", edgecolor="#bdbdbd", lw=0.4)
    gt.plot(ax=axn, facecolor="#d9d9d9", edgecolor="#6b6b6b", lw=0.6)
    centro = zonas_est.to_crs("EPSG:4326").union_all().centroid
    axn.plot(centro.x, centro.y, marker="s", ms=5, mfc="#d73027", mec=TINTA, mew=0.6)
    axn.text(centro.x + 0.15, centro.y - 0.05, "Ciudad de\nGuatemala", fontsize=5.5, va="top", color=TINTA)
    encuadre(axn, gt.total_bounds, 89.4, 50, margen=0.15)
    axn.set_facecolor("#f7fbff")
    marco(axn, "Ubicación nacional")

    # Ubicación municipal
    axm = fig.add_axes(rect(180, 80, 89.4, 50))
    zonas_mun.plot(ax=axm, facecolor="#e0e0e0", edgecolor="#9e9e9e", lw=0.2)
    zonas_est.plot(ax=axm, facecolor="#d73027", edgecolor=TINTA, lw=0.4, alpha=0.85)
    encuadre(axm, zonas_mun.total_bounds, 89.4, 50)
    marco(axm, "Ubicación municipal – municipio de Guatemala")

    # Leyenda
    axl = fig.add_axes(rect(180, 134, 89.4, 38))
    marco(axl)
    axl.text(0.04, 0.93, "LEYENDA", transform=axl.transAxes, fontsize=7.5, fontweight="bold", va="top")
    leyenda(axl)

    # Cuadro técnico
    axc = fig.add_axes(rect(180, 175, 89.4, 31))
    marco(axc)
    lineas = ["Sistema de referencia: UTM 15N (EPSG:32615) · WGS84"] + tecnico_extra + [
        f"Período: {periodo_txt}", f"Elaboración: {AUTOR}"]
    axc.text(0.04, 0.92, "\n".join(lineas), transform=axc.transAxes, fontsize=6, va="top",
             linespacing=1.55, color=TINTA)

    fig.savefig(archivo.with_suffix(".pdf"))
    fig.savefig(archivo.with_suffix(".png"), dpi=300)
    plt.close(fig)


def entradas_leyenda(axl, items, y0=0.78, dy=0.105):
    for i, (tipo, color, texto) in enumerate(items):
        y = y0 - i * dy
        if tipo == "parche":
            axl.add_patch(mpatches.Rectangle((0.05, y - 0.035), 0.09, 0.07, transform=axl.transAxes,
                                             facecolor=color, edgecolor="#6b6b6b", lw=0.4))
        elif tipo == "linea":
            axl.add_line(Line2D([0.05, 0.14], [y, y], transform=axl.transAxes, color=color, lw=1.0))
        axl.text(0.18, y, texto, transform=axl.transAxes, va="center", fontsize=6.5)


def barra_color(axl, cmap, norm, etiqueta, ticks=None, y=0.30):
    cax = axl.inset_axes([0.06, y, 0.86, 0.09])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = plt.colorbar(sm, cax=cax, orientation="horizontal", ticks=ticks)
    cb.outline.set_linewidth(0.4)
    cb.ax.tick_params(labelsize=6, length=2)
    axl.text(0.06, y + 0.13, etiqueta, transform=axl.transAxes, fontsize=6.5)


# ============================ MAPAS ============================
SALIDA.mkdir(parents=True, exist_ok=True)
generados = []


def quiere(nombre):
    return SOLO is None or nombre in SOLO


# --- 1. ICU por intensidad ---
COL_ICU = ["#2166ac", "#f7f7f7", "#fdae61", "#d73027"]
cmap_icu = ListedColormap(COL_ICU)
norm_icu = BoundaryNorm([0.5, 1.5, 2.5, 3.5, 4.5], 4)
for clave, (largo, tif, _, ptxt) in PERIODOS.items():
    nombre = f"ICU_{clave}"
    if not quiere(nombre):
        continue
    a, ext = leer_raster(BASE / "ICU_Clases" / "Icu_Clases" / f"ICU_clases_{tif}.tif")
    a[(a < 1) | (a > 4)] = np.nan

    def dibujar(ax, a=a, ext=ext):
        ax.imshow(a, extent=ext, cmap=cmap_icu, norm=norm_icu, interpolation="nearest", zorder=5)

    def leyenda(axl):
        axl.text(0.05, 0.80, "Clases de intensidad térmica", transform=axl.transAxes, fontsize=6.5,
                 fontstyle="italic")
        entradas_leyenda(axl, [("parche", COL_ICU[0], "Fría  (< μ − σ)"),
                               ("parche", COL_ICU[1], "Templada  (μ − σ a μ + σ)"),
                               ("parche", COL_ICU[2], "Cálida / ICU  (> μ + σ)"),
                               ("parche", COL_ICU[3], "Muy cálida  (> μ + 2σ)"),
                               ("linea", TINTA, "Límite de zonas 1, 4, 9 y 10")], y0=0.66, dy=0.125)

    lamina(f"DISTRIBUCIÓN TÉRMICA SUPERFICIAL POR INTENSIDAD – {largo.upper()}", dibujar, leyenda, ptxt,
           ["Fuente: Landsat 8/9 C2 (USGS), Google Earth Engine", f"Método: {METODO_LST}",
            "Clases: μ ± σ de la LST de cada período"], SALIDA / f"Mapa_{nombre}")
    generados.append(nombre)

# --- 2. LST estacional (escala común) ---
lst = {k: leer_raster(BASE / "LST_SplitWindow" / f"LST_SplitWindow_{v[1]}.tif") for k, v in PERIODOS.items()}
vmin = np.floor(min(np.nanpercentile(a, 1) for a, _ in lst.values()))
vmax = np.ceil(max(np.nanpercentile(a, 99) for a, _ in lst.values()))
cmap_lst = plt.get_cmap("YlOrRd")
for clave, (largo, tif, _, ptxt) in PERIODOS.items():
    nombre = f"LST_{clave}"
    if not quiere(nombre):
        continue
    a, ext = lst[clave]

    def dibujar(ax, a=a, ext=ext):
        ax.imshow(a, extent=ext, cmap=cmap_lst, vmin=vmin, vmax=vmax, interpolation="nearest", zorder=5)

    def leyenda(axl, a=a):
        barra_color(axl, cmap_lst, plt.Normalize(vmin, vmax), "Temperatura superficial (°C)", y=0.52)
        axl.text(0.06, 0.36, f"Media del área: {np.nanmean(a):.1f} °C\n"
                             f"Mín.–máx.: {np.nanmin(a):.1f}–{np.nanmax(a):.1f} °C",
                 transform=axl.transAxes, fontsize=6.2, va="top")
        entradas_leyenda(axl, [("linea", TINTA, "Límite de zonas 1, 4, 9 y 10")], y0=0.12)

    lamina(f"TEMPERATURA SUPERFICIAL TERRESTRE (LST) – {largo.upper()}", dibujar, leyenda, ptxt,
           ["Fuente: Landsat 8/9 C2 (USGS), Google Earth Engine", f"Método: {METODO_LST}",
            "Compuesto: mediana estacional · escala común 4 períodos"], SALIDA / f"Mapa_{nombre}")
    generados.append(nombre)

# --- 3. NDVI estacional ---
cmap_ndvi = plt.get_cmap("Greens")
for clave, (largo, tif, tif_s2, ptxt) in PERIODOS.items():
    nombre = f"NDVI_{clave}"
    if not quiere(nombre):
        continue
    ruta = BASE / "Variables- NDVI NDBI Atmos" / f"S2_NDVI_NDBI_{tif_s2}.tif"
    if not ruta.exists():
        print(f"Falta {ruta.name}; se omite {nombre}")
        continue
    a, ext = leer_raster(ruta, banda=1)

    def dibujar(ax, a=a, ext=ext):
        ax.imshow(a, extent=ext, cmap=cmap_ndvi, vmin=0, vmax=0.8, interpolation="nearest", zorder=5)

    def leyenda(axl, a=a):
        barra_color(axl, cmap_ndvi, plt.Normalize(0, 0.8), "NDVI (B8 − B4) / (B8 + B4)",
                    ticks=[0, 0.2, 0.4, 0.6, 0.8], y=0.52)
        axl.text(0.06, 0.36, f"NDVI medio del área: {np.nanmean(a):.2f}", transform=axl.transAxes,
                 fontsize=6.2, va="top")
        entradas_leyenda(axl, [("linea", TINTA, "Límite de zonas 1, 4, 9 y 10")], y0=0.12)

    lamina(f"ÍNDICE DE VEGETACIÓN DE DIFERENCIA NORMALIZADA (NDVI) – {largo.upper()}", dibujar, leyenda,
           ptxt, ["Fuente: Sentinel-2 L2A (ESA/Copernicus), Google Earth Engine",
                  "Compuesto: mediana estacional · 10 m"], SALIDA / f"Mapa_{nombre}")
    generados.append(nombre)

# --- 4. Relación LST–NDVI: coeficiente local MGWR ---
ruta_mgwr = next((p for p in [BASE / "OLS_GWR" / "malla_modelo_local.gpkg",
                   BASE / "GEE_ICU" / "OLS_GWR" / "malla_modelo_local.gpkg"] if p.exists()),
                  BASE / "OLS_GWR" / "malla_modelo_local.gpkg")
if ruta_mgwr.exists():
    malla = gpd.read_file(ruta_mgwr).to_crs(CRS)
    cols = [f"NDVI_degC_{k}" for k in PERIODOS]
    lim = float(np.nanpercentile(np.abs(malla[cols].values), 99))
    cmap_c = plt.get_cmap("RdBu_r")
    norm_c = TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
    for clave, (largo, tif, _, ptxt) in PERIODOS.items():
        nombre = f"MGWR_NDVI_{clave}"
        if not quiere(nombre):
            continue
        col = f"NDVI_degC_{clave}"
        sig = malla[f"sig_NDVI_{clave}"] == 1

        def dibujar(ax, col=col, sig=sig):
            malla[~sig].plot(ax=ax, facecolor="#d9d9d9", edgecolor="none", zorder=5)
            malla[sig].plot(ax=ax, column=col, cmap=cmap_c, norm=norm_c, edgecolor="none", zorder=5)

        def leyenda(axl, col=col, sig=sig):
            barra_color(axl, cmap_c, norm_c, "Cambio de LST por unidad de NDVI (°C)", y=0.55)
            axl.text(0.06, 0.40, f"Mediana: {np.median(malla[col]):.1f} °C · "
                                 f"celdas significativas: {100 * sig.mean():.0f} %",
                     transform=axl.transAxes, fontsize=6.2, va="top")
            entradas_leyenda(axl, [("parche", "#d9d9d9", "No significativo (α = 0.05, t corregido)"),
                                   ("linea", TINTA, "Límite de zonas 1, 4, 9 y 10")], y0=0.24, dy=0.12)

        lamina(f"RELACIÓN LOCAL LST–NDVI (MGWR) – {largo.upper()}", dibujar, leyenda, ptxt,
               ["Modelo: MGWR LST ~ NDVI + NDBI · malla 100 m", "Kernel bisquare adaptativo · AICc",
                "Fuente: Landsat 8/9 C2 y Sentinel-2 L2A"], SALIDA / f"Mapa_{nombre}")
        generados.append(nombre)

print("Mapas generados:", ", ".join(generados))
print(f"Carpeta: {SALIDA}")
