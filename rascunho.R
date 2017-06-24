

#ff = fB*bin_pos
#display(ff/max(ff))
#edgeff = edge.detect(ff,thresh1=1, thresh2=15, noise="gaussian", noise.s=3,
#                     method="Canny")
#display(edgeff)
#fulledge = fillHull(edgeff)
#display(fulledge)

# Vamos implememtar um seedgrow ? pq?

# Reconstroi usando so a aproximada da DWT
#wave = dwt_matrix(normalizada)
#display(wave/max(wave))
#wave = extrai1Q(wave)
#display(wave/max(wave))
#wave2 = idwt_matrix(wave)
#display(wave2/max(wave2))
#wave2=normalizada

#maximoLocal = diff(sign(diff(h$y) ) )==-2
#picos = h$x[ (which(significativo & maximoLocal) )+1 ] 
#p1 = picos[1]
#p2 = picos[2]
#### Encontrar intervalo de valores significativos #