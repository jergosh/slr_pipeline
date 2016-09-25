data=read.csv2(file="fichierdepartements.csv")
x=as.numeric(as.vector(data$Latitude))
y=as.numeric(as.vector(data$Longitude))

vec=read.table("data.txt")
v=as.vector(vec$V1)

dist=function(x1,y1,x2,y2) {(x1-x2)^2+(y1-y2)^2}


n=nrow(vec) # Correct?
print(n)
maxconc=0
maxcenter=c(0,0)
maxrad=0
maxvec=0
moy=mean(v)


for (i in 1:n) {
  for (j in 1:n) {
    SV=0
    m=0
    rad=dist(x[i],x[j],y[i],y[j])
    vec=0
    for (k in 1:n) {
      if (dist(x[i],x[k],y[i],y[k])<=rad) {
        m=m+1
        SV=SV+v[k]
        vec=c(vec,k)
      }
    }
    conc=0
    if (m<n) conc=(SV-m*moy)/sqrt(m*(n-m))    
    if (conc>maxconc) {
      maxconc=conc
      maxcenter=c(x[i],y[i])
      maxrad=rad
      maxvec=vec
    }
    
  }
}





Elections:
  > maxconc
[1] 3.501
> maxvec
[1]  0  2  3  9 11 12 14 15 16 17 18 19 21 22 23 26 27 28 30 31 32 34 35 36 39
[26] 40 43 44 45 46 48 49 52 55 57 58 59 60 61 62 63 64 71 74 75 76 77 78 79 80
[51] 81 84 85 86 88 90 91 92 93 94
p-val=0.001

Revenu moyen:
  > maxconc
[1] 2410.952
> maxvec
[1]  0 74 77 91
pval=0.001


n=94
minconc=0
mincenter=c(0,0)
minrad=0
minvec=0

for (i in 1:n) {
  for (j in 1:n) {
    SV=0
    m=0
    rad=dist(x[i],x[j],y[i],y[j])
    vec=0
    for (k in 1:n) {
      if (dist(x[i],x[k],y[i],y[k])<=rad) {
        m=m+1
        SV=SV+v[k]
        vec=c(vec,k)
      }
    }
    conc=0
    if (m<n) conc=(SV-m*moy)/sqrt(m*(n-m))    
    if (conc<minconc) {
      minconc=conc
      mincenter=c(x[i],y[i])
      minrad=rad
      minvec=vec
    }
    
  }
}



Elections:
  > minconc
[1] -3.092210
> minvec
[1]  0  1  3  4  5  6  7 10 13 20 24 25 29 33 37 38 41 42 47 50 51 53 54 56 57
[26] 62 66 67 68 69 70 72 73 82 83 87 88 89

Revenu moyen:
  > minconc
[1] -1357.523
> minvec
[1]  0  3  7  9 11 12 15 16 17 18 19 22 23 25 29 30 31 32 33 35 39 41 42 45 46
[26] 47 57 62 63 64 65 78 80 81 83 85 86
pval=0.126


n=94
maxconcobs=2410.952
minconcobs=-1357.23
nsim=999
pmin=0
pmax=0
for (sim in 1:nsim) {
  print(sim)
  vsim=v[sample(n)]
  maxconc=0
  minconc=0
  
  for (i in 1:n) {
    for (j in 1:n) {
      SV=0
      m=0
      rad=dist(x[i],x[j],y[i],y[j])
      for (k in 1:n) {
        if (dist(x[i],x[k],y[i],y[k])<=rad) {
          m=m+1
          SV=SV+vsim[k]
        }
      }
      conc=0
      if (m<n) conc=(SV-m*moy)/sqrt(m*(n-m))        
      if (conc>maxconc) maxconc=conc
      if (conc<minconc) minconc=conc
      
    }
  }
  
  if (maxconc>maxconcobs) pmax=pmax+1
  if (minconc<minconcobs) pmin=pmin+1
  
  
}






----------------------------------------------------
  
  
  
 data=read.csv2(file="fichierdepartements.csv")
x=as.numeric(as.vector(data$Latitude))
y=as.numeric(as.vector(data$Longitude))
vec=read.table("datacas.txt")
veccas=as.vector(vec$V1)
vec=read.table("datacontroles.txt")
veccontroles=as.vector(vec$V1)


dist=function(lat1,lat2,lon1,lon2) {(lat1-lat2)^2+(lon1-lon2)^2}


n=94
maxconc=0
maxcenter=c(0,0)
maxrad=0
maxvec=0
totalcontroles=sum(veccontroles)
moy=sum(veccas)/totalcontroles


for (i in 1:n) {
  for (j in 1:n) {
    Scas=0
    Scontroles=0
    rad=dist(x[i],x[j],y[i],y[j])
    vec=0
    for (k in 1:n) {
      if (dist(x[i],x[k],y[i],y[k])<=rad) {
        Scas=Scas+veccas[k]
        Scontroles=Scontroles+veccontroles[k]
        vec=c(vec,k)
      }
    }
    conc=0
    if (Scontroles<totalcontroles) conc=(Scas-Scontroles*moy)/sqrt(Scontroles*(totalcontroles-Scontroles))    
    if (conc>maxconc) {
      maxconc=conc
      maxcenter=c(x[i],y[i])
      maxrad=rad
      maxvec=vec
    }
    
  }
}


Elections:
  > maxconc
[1] 0.03878143
> maxvec
[1]  0  1  2  3  7  8  9 10 11 12 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28
[26] 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 51 52 53
[51] 54 55 56 57 58 59 60 61 62 63 64 65 68 69 70 71 72 74 75 76 77 78 79 80 81
[76] 84 85 86 87 88 90 91 92 93 94
p-val=0.001

Revenu moyen:
  > maxconc
[1] 3714.678
> maxvec
[1]  0 74 77 91
p-val=0.002


n=94
minconc=0
mincenter=c(0,0)
minrad=0
minvec=0
totalcontroles=sum(veccontroles)
moy=sum(veccas)/totalcontroles


for (i in 1:n) {
  for (j in 1:n) {
    Scas=0
    Scontroles=0
    rad=dist(x[i],x[j],y[i],y[j])
    vec=0
    for (k in 1:n) {
      if (dist(x[i],x[k],y[i],y[k])<=rad) {
        Scas=Scas+veccas[k]
        Scontroles=Scontroles+veccontroles[k]
        vec=c(vec,k)
      }
    }
    conc=0
    if (Scontroles<totalcontroles) conc=(Scas-Scontroles*moy)/sqrt(Scontroles*(totalcontroles-Scontroles))    
    if (conc<minconc) {
      minconc=conc
      mincenter=c(x[i],y[i])
      minrad=rad
      minvec=vec
    }
    
  }
}


Elections:
  > minconc
[1] -0.03406365
> minvec
[1]  0  1  3  4  5  6  7 10 13 20 24 25 29 33 37 38 41 42 47 50 51 53 54 56 57
[26] 62 66 67 68 69 70 72 73 82 83 87 88 89
p-val=0.001

Revenu moyen:
  > minconc
[1] -1718.523
> minvec
[1]  0  3  7  9 11 12 13 14 15 16 17 18 19 21 22 23 25 26 27 28 29 30 31 32 33
[26] 34 35 36 39 40 41 42 43 44 45 46 47 48 49 52 55 57 60 62 63 64 65 68 70 71
[51] 78 80 81 83 84 85 86 88
p-val=0.492


n=94
maxconcobs=3714.678
minconcobs=-1718.523
nsim=999
pmin=0
pmax=0
for (sim in 1:nsim) {
  print(sim)
  ordsim=sample(n)    
  maxconc=0
  minconc=0
  
  for (i in 1:n) {
    for (j in 1:n) {
      Scas=0
      Scontroles=0
      rad=dist(x[i],x[j],y[i],y[j])
      for (k in 1:n) {
        if (dist(x[i],x[k],y[i],y[k])<=rad) {
          Scas=Scas+veccas[ordsim[k]]
          Scontroles=Scontroles+veccontroles[ordsim[k]]
        }
      }
      conc=0
      if (Scontroles<totalcontroles) conc=(Scas-Scontroles*moy)/sqrt(Scontroles*(totalcontroles-Scontroles))            
      if (conc>maxconc) maxconc=conc
      if (conc<minconc) minconc=conc
      
    }
  }
  
  if (maxconc>maxconcobs) pmax=pmax+1
  if (minconc<minconcobs) pmin=pmin+1
  
  
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  Cancer:
  
library(spatstat)
data(chorley)
x=chorley$x
y=chorley$y
v=as.numeric(chorley$marks=="larynx")

dist=function(lat1,lat2,lon1,lon2) {(lat1-lat2)^2+(lon1-lon2)^2}


n=1036
maxconc=0
maxcenter=c(0,0)
maxrad=0
maxvec=0
moy=mean(v)


for (i in 1:58) {
  for (j in 1:58) {
    SV=0
    m=0
    rad=dist(x[i],x[j],y[i],y[j])
    vec=0
    for (k in 1:n) {m
                    if (dist(x[i],x[k],y[i],y[k])<=rad) {
                      m=m+1
                      SV=SV+v[k]
                      vec=c(vec,k)
                    }
    }
    conc=0
    if (m<n) conc=(SV-m*moy)/sqrt(m*(n-m))    
    if (conc>maxconc) {
      maxconc=conc
      maxcenter=c(x[i],y[i])
      maxrad=rad
      maxvec=vec
    }
    
  }
}

> maxconc
[1] 0.05181287
> maxvec
[1]   0  55  56  57  58 278



pval=1
nsim=999

for (sim in 1:nsim) {
  print (sim)
  perm=sample(n)
  xsim=x[perm]
  ysim=y[perm]
  maxconcsim=0
  
  for (i in 1:58) { # This seems to be an optimisation, to only iterate over non-zero marks
    for (j in 1:58) {
      SV=0
      m=0
      rad=dist(xsim[i],xsim[j],ysim[i],ysim[j])
      vec=0
      for (k in 1:n) {
        m
        if (dist(xsim[i],xsim[k],ysim[i],ysim[k])<=rad) {
          m=m+1
          SV=SV+v[k]
        }
      }
      conc=0
      if (m<n)
        conc=(SV-m*moy)/sqrt(m*(n-m))    
      if (conc>maxconcsim) {
        maxconcsim=conc
      }
    }
  }

  print(maxconcsim)
  if (maxconcsim>=maxconc)
    pval=pval+1
}

pval=pval/(nsim+1)    

> pval
[1] 0.018

dont 10 ex-aequo => pval=0.013
