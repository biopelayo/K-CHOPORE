# K-CHOPORE v2.0 — Informe General

### Qué hay de nuevo y por qué importa

**Autor:** Pelayo Gonzalez de Lena Rodriguez
**Fecha:** Febrero 2026
**Proyecto:** Análisis de RNA directo por nanoporos en *Arabidopsis thaliana*

---

## En una frase

K-CHOPORE v2 convierte datos crudos de secuenciación de nanoporos en resultados biológicos listos para publicar — de forma automática, reproducible, y sin necesidad de tocar la línea de comandos más de una vez.

---

## 1. ¿Qué es K-CHOPORE?

Imagina que tienes una planta (*Arabidopsis*) y quieres saber **qué genes están encendidos o apagados** cuando la sometes a un estrés (en nuestro caso, una droga llamada Antimicina A que bloquea la respiración de las mitocondrias).

Para eso, extraes el RNA de la planta (las "fotocopias" de los genes activos) y lo pasas por un secuenciador de nanoporos (Oxford Nanopore MinION) — un aparato del tamaño de un USB que lee cada molécula de RNA individualmente.

El problema es que **del secuenciador salen millones de datos crudos** que hay que procesar con muchos programas distintos, en un orden concreto, con los parámetros correctos. Eso es lo que hace K-CHOPORE: **automatiza todo el análisis**, desde los datos crudos hasta las tablas de genes y los gráficos finales.

---

## 2. ¿Qué había antes (v1)?

La versión original era un **prototipo funcional sobre el papel** pero con problemas importantes en la práctica:

- **No se podía ejecutar sin modificarlo a mano** — había errores en las instrucciones que recibían los programas internos
- **Faltaban programas** — algunas herramientas clave no estaban incluidas en el paquete
- **No había forma de probar que funcionaba** — no existían tests automáticos
- **El análisis estadístico era básico** — solo comparaba dos grupos simples, no podía manejar un diseño experimental complejo
- **No tenía documentación práctica** — un usuario nuevo no sabía por dónde empezar

En resumen: era un borrador prometedor que necesitaba mucho trabajo para funcionar con datos reales.

---

## 3. ¿Qué tiene de nuevo K-CHOPORE v2?

### 3.1 Ahora funciona de verdad ✅

Suena obvio, pero es el cambio más importante. Se corrigieron **6 errores críticos** que impedían que el pipeline funcionara con datos reales:

| Problema | Consecuencia | Solución |
|----------|-------------|----------|
| El programa de alineamiento (minimap2) usaba el modo genérico para nanoporos | Las lecturas de RNA se alineaban en la dirección equivocada y sin reconocer los empalmes de los genes | Se activó el modo splice-aware con orientación directa RNA |
| Los nombres de los cromosomas no coincidían entre el genoma de referencia y la anotación de genes | El análisis de isoformas producía cero resultados | Se añadió traducción automática de nombres |
| El programa de isoformas (FLAIR) fallaba silenciosamente con guiones bajos en los nombres | Los resultados se asignaban a las muestras equivocadas | Se detectó y corrigió el formato de nombres |
| El análisis estadístico no reconocía las columnas de datos | DESeq2 no podía emparejar muestras con conteos | Se generan automáticamente las correspondencias correctas |

**Resultado:** El pipeline procesó con éxito **10 muestras × 78+ pasos** sin intervención manual.

### 3.2 Todo dentro de una caja Docker 📦

Docker es como una "cápsula" que contiene todo lo necesario para ejecutar el pipeline: los 15+ programas bioinformáticos, las versiones exactas de cada librería, y el sistema operativo. Ventajas:

- **Reproducibilidad total** — Cualquier persona, en cualquier ordenador, obtiene exactamente los mismos resultados
- **Sin conflictos de instalación** — No hace falta instalar nada más que Docker
- **Funciona en Windows, Mac y Linux** — Probado en los tres sistemas

La imagen Docker pesa 22.7 GB porque incluye TODO: desde el alineador hasta el programa de estadística, pasando por los detectores de modificaciones del RNA.

### 3.3 Diseño experimental de verdad 🔬

La versión original solo podía comparar "Grupo A vs Grupo B". La v2 maneja un **diseño factorial 2×2 completo**:

```
                    Control     Antimicina A
    WT (normal)     3 réplicas  3 réplicas
    anac017-1       3 réplicas  1 réplica
    (mutante)
```

Esto permite responder tres preguntas simultáneamente:
1. **Efecto del genotipo:** ¿Qué genes cambian por ser mutante?
2. **Efecto del tratamiento:** ¿Qué genes cambian con la droga?
3. **Interacción:** ¿Responde el mutante de forma distinta a la droga que el normal?

Además, el programa detecta automáticamente si algún grupo tiene pocas réplicas y ajusta el modelo estadístico para no dar resultados falsos (en nuestro caso, usa un modelo aditivo porque el grupo mutante+droga solo tiene 1 réplica).

### 3.4 Análisis de isoformas (no solo genes) 🧬

La mayoría de pipelines de RNA-seq cuentan cuántas lecturas tiene cada **gen**. K-CHOPORE va un paso más allá: cuenta **isoformas** — las distintas versiones (empalmes alternativos) de cada gen.

Esto es posible porque la secuenciación de nanoporos lee cada molécula de RNA de principio a fin (lecturas de ~1000 nucleótidos de media), mientras que las tecnologías convencionales cortan el RNA en trocitos de 150 nt y luego intentan reconstruir el puzzle.

**Resultado:** Se cuantificaron **20.958 isoformas** distintas en nuestras 10 muestras.

### 3.5 Detección de modificaciones del RNA 🏷️

El RNA no es solo una secuencia de letras (A, U, G, C) — lleva "marcas" químicas que regulan su función. La secuenciación directa de nanoporos puede detectar estas marcas porque alteran la señal eléctrica cuando el RNA pasa por el poro.

K-CHOPORE integra dos herramientas de detección:
- **ELIGOS2** — detecta modificaciones analizando los errores de lectura (no requiere datos de señal)
- **m6Anet** — detecta metilaciones m6A a partir de la señal eléctrica original

> *Nota: en esta ejecución, ELIGOS2 falló por un problema de compatibilidad software conocido, y m6Anet no se ejecutó porque los archivos de señal (644 GB) no cabían en el disco del servidor. Ambos quedan como funcionalidad disponible para futuras ejecuciones.*

### 3.6 Informes de calidad automáticos 📊

K-CHOPORE genera automáticamente:

- **NanoPlot** — Un informe por cada muestra con histogramas de longitud, calidad, y gráficos de dispersión
- **NanoComp** — Comparación visual entre todas las muestras (para detectar muestras problemáticas)
- **pycoQC** — Informe interactivo de la calidad del secuenciador (rendimiento del flow cell, reads por canal, etc.)
- **MultiQC** — Un informe resumen que junta toda la información de calidad en una sola página web

Todo esto se genera sin intervención manual.

### 3.7 Tests automáticos 🧪

Se creó un conjunto de **51 tests automáticos** que verifican:
- Que la configuración del pipeline sea correcta
- Que minimap2 use los parámetros adecuados para RNA directo
- Que todas las reglas del pipeline estén conectadas correctamente
- Que los datos de prueba produzcan los resultados esperados

Esto permite que cualquier modificación futura se pueda validar al instante.

### 3.8 Documentación completa 📖

- **Guía paso a paso** — Tutorial de 548 líneas para un usuario que nunca ha usado el pipeline
- **Sección de Métodos** — Texto listo para copiar-pegar en un artículo científico
- **README actualizado** — Con ejemplos reales, troubleshooting, y estructura del proyecto

---

## 4. Herramientas integradas

| Herramienta | Función | Estado |
|-------------|---------|--------|
| **NanoFilt** | Filtrar lecturas de baja calidad | ✅ Funcional |
| **NanoPlot** | Gráficos de calidad por muestra | ✅ Funcional |
| **NanoComp** | Comparación entre muestras | ✅ Funcional |
| **pycoQC** | Calidad del secuenciador | ✅ Funcional |
| **minimap2** | Alineamiento splice-aware | ✅ Funcional (corregido) |
| **samtools** | Estadísticas de alineamiento | ✅ Funcional |
| **FLAIR** | Detección y cuantificación de isoformas | ✅ Funcional (corregido) |
| **ELIGOS2** | Modificaciones del RNA (error-based) | ⚠️ Fallo CMH test |
| **m6Anet** | Detección de m6A (signal-level) | ⏸️ Requiere FAST5 en disco |
| **Nanopolish** | Análisis de señal (para m6Anet) | ✅ Compilado correctamente |
| **DESeq2** | Expresión diferencial factorial | ✅ Funcional (reescrito) |
| **MultiQC** | Informe resumen de calidad | ✅ Funcional |
| **StringTie2** | Ensamblaje de transcritos (disponible) | ⏸️ No usado en este análisis |

---

## 5. Resultados obtenidos con K-CHOPORE v2

### 5.1 Datos procesados

| Métrica | Valor |
|---------|-------|
| Muestras analizadas | 10 |
| Lecturas crudas totales | 12.763.204 |
| Lecturas tras filtro de calidad | 12.603.968 (98,8%) |
| Longitud media de lectura | 997 nucleótidos |
| Tasa media de alineamiento | 90,5% |
| Isoformas cuantificadas | 20.958 |

### 5.2 Genes diferencialmente expresados

**Por genotipo** (mutante vs normal):
- 435 isoformas con cambio significativo
- 303 activadas en el mutante, 132 reprimidas
- Gen más afectado: un transportador de magnesio mitocondrial (AT1G53480) — 18 veces menos expresado en el mutante

**Por tratamiento** (Antimicina A vs control):
- 266 isoformas con cambio significativo
- El 92% están activadas (la droga enciende genes de defensa)
- **AOX1A** (Oxidasa Alternativa 1A) aparece entre los más activados — esto es importantísimo porque AOX1A es el marcador clásico de la respuesta mitocondrial, lo que **valida que todo el experimento y el análisis funcionan correctamente**

### 5.3 Gráficos generados

- **PCA** — Muestra que el 44% de la variabilidad se debe al genotipo y el 28% al tratamiento
- **Volcano plots** — Visualizan qué genes cambian más y cuáles son estadísticamente significativos
- **MA plots** — Muestran la relación entre la cantidad de expresión y el cambio observado
- **Violín plots** — Comparan distribuciones de calidad y longitud entre muestras

---

## 6. ¿Qué queda por hacer?

| Tarea | Dificultad | Impacto |
|-------|-----------|---------|
| Resolver el fallo de ELIGOS2 (actualizar rpy2/R) | Media | Alto — permitiría detectar modificaciones del RNA |
| Ejecutar m6Anet (necesita espacio en disco para 644 GB de FAST5) | Baja (logística) | Alto — detección de metilación m6A |
| Completar las 2 réplicas faltantes del grupo mutante×AA (requiere GPU para basecalling) | Media | Medio — permitiría estimar la interacción genotipo×tratamiento |
| Integrar análisis funcional (GO enrichment) | Baja | Medio — contexto biológico de los genes diferenciales |

---

## 7. Resumen visual

```
ANTES (v1)                           DESPUÉS (v2)
─────────────────────────────────    ─────────────────────────────────
❌ Errores críticos que impedían     ✅ 10 muestras procesadas con
   la ejecución                         éxito (78+ pasos automáticos)

❌ No se podía reproducir            ✅ Docker: mismos resultados
                                        en cualquier ordenador

❌ Comparación simple A vs B         ✅ Diseño factorial 2×2 con
                                        modelo adaptativo

❌ Solo contaba genes                ✅ 20.958 isoformas cuantificadas

❌ Sin tests                         ✅ 51 tests automáticos

❌ Sin documentación práctica        ✅ Guías, métodos, troubleshooting

❌ ELIGOS2 no instalado              ✅ ELIGOS2 instalado y parcheado
                                        (fallo CMH pendiente)

❌ Nanopolish no compilaba           ✅ Compilación correcta (serial)

❌ Sin scripts de despliegue         ✅ Despliegue con un solo comando
```

---

## 8. Cómo ejecutar K-CHOPORE v2

Para quien quiera replicar el análisis:

```bash
# 1. Clonar el repositorio
git clone https://github.com/biopelayo/K-CHOPORE.git

# 2. Construir la imagen Docker (una sola vez, ~45 min)
docker build -t k-chopore:latest .

# 3. Configurar las muestras en config/config.yml

# 4. Ejecutar
docker run --rm -v $(pwd):/workspace -w /workspace k-chopore:latest \
    snakemake --cores 40 --keep-going --latency-wait 60
```

Los resultados aparecen en la carpeta `results/` con todos los informes, gráficos y tablas.

---

*K-CHOPORE v2.0 — Pipeline de análisis de RNA directo por nanoporos*
*Universidad de Oviedo | FINBA | Febrero 2026*
