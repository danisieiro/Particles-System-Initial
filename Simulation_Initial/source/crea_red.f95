!   PROGRAMA SIMULACIÓN CREA_RED
!
!   ESTE PROGRAMA CREA UNA CONFIGURACIÓN INICIAL DE UN SISTEMA DE N PARTÍCULAS EN UNA CAJA DE LADO L.
!   ESTAS PARTÍCULAS SE COLOCAN EN FORMA DE RED FCC (FACE CENTERED CUBIC) CON PARAMETRO DE RED L/k (NUESTRO CASO K=5)
!   DE FORMA QUE TENDREMOS N/4 "CAJAS" QUE CONTIENEN 4 PARTÍCULAS CADA UNA.
!
!   MEDIANTE UNA FUNCION QUE DEVUELVE NÚMEROS ALEATORIOS ENTRE 0 Y 1, INTRODUCIMOS UNA MODIFICACIÓN
!   SOBRE CADA PARTICULA DE HASTA UN 1% DEL LADO L. INTRODUCIMOS UNA VELOCIDAD A CADA PARTÍCULA
!   ENTRE 0 Y 1, PERMITIÉNDONOS CALCULAR LA ecin.
!
!   SE CALCULA LA ENERGÍA POTENCIAL ASUMIENDO UN POTENCIAL DE LENNARD-JONES ENTRE PARES DE PARTÍCULAS
!   OBTENEMOS ENERGIA POTENCIAL, ACELERACIONES Y DERIVADAS DE FI RESPECTO DE V
!   
!   FIJAMOS ENERGÍA TOTAL (MI CASO etot = -615) Y LE RESTAMOS LA POTENCIAL, OBTENIENDO ECINQ,
!   MULTIPLICAMOS VELOCIDADES POR EL FACTOR (ECINQ/ecin)**(-1/2).
!
!   GRABAMOS PARÁMETROS, POSICIONES, VELOCIDADES Y ACELERACIONES
!

program crea_red

    use def_precision


    implicit none

!   np = indice de la partícula (va de 1 a npmax/(k-1))
!   pl = lado de la caja
!   vol = volumen, pl*pl*pl
!   numk = cajitas o porciones en las que dividimos la caja de lado pl
!   k3 = numero de cajitas totales que forman el volumen total
!   npmax = numero de particulas maxima (parametro), en este caso 500
!   pa = lado de cada una de las cajitas pequeñas
!   pma = mitad de pma
!   rc = radio de corte, a partir del cual no se consideran interacciones entre pares (mitad del lado de la caja)
!   rc2 = radio de corte al cuadrado
!   dens = densidad de particulas
!
!   rx,rx,rz = array de posiciones en las coordenadas x,y,z
!   vx,vy,vz = array de velocidades en las direcciones x,y,z
!   ax,ay,az = array de aceleraciones en las direcciones x,y,z
!
!   px,py,pz = suma de velocidades vz,vy,vz
!
!   epot = energía potencial total
!   ecin = energia cinetica total
!   etot = suma de energias cinérica y potencial
!   etotq = energía total que queremos fijar
!   ecinq = energía cinética fijada luego de fijar la energia total (epot no cambia)
!
!   indices de bucles: i, j, k, l
!
!   dfiv = derivada del potencial con respecto al volumen
!   d2fiv = derivada segunda del potencial con respecto al volumen
!
!   fname = archivo de parametros
!   gname = archivo de posiciones, velocidades y aceleraciones en binario



!   DECLARAMOS LAS VARIABLES
                
    integer (kind = entero) :: np,pl,vol,numk,K3
    integer (kind=entero), parameter :: npmax=500
    real (kind = doblep), parameter :: pi=3.141592653589
    real (kind = doblep), dimension(npmax) :: rx,ry,rz,vx,vy,vz,ax,ay,az
    real (kind = doblep) :: epot,dfiv,d2fiv
    integer (kind = entero) :: i,j,k,l
    real (kind = doblep) :: pma,iflag,px,py,pz,mi_random,rc,rc2,pli,xnp,dens,pa,ecinq,ecin,etotq,factor
    character(len=50) :: fname, gname
    
!   INTRODUCIR LOS VALORES NECESARIOS
    numk = 5
    k3 = numk*numk*numk
    xnp=dble(npmax)
    pl = 10
    pli=1.d00/pl
    vol=pl*pl*pl
    dens=xnp/dble(vol)
    rc = dble(pl/2.d00)
    rc2 = rc * rc
    pa=dble(pl/numk)
    pma=dble(pa)/2.d00



!   PONEMOS EL INDICE DE PARTICULAS A 0
    np=0


!   COLOCO LAS PARTÍCULAS EN SUS POSICIONES, YENDO 'CAJA A CAJA' (5X5X5 CAJAS CON 4 PARTÍCULAS EN CADA UNA)

    do i=0,numk-1
        do j=0, numk-1
            do k=0, numk-1
!               SIGUIENTE PARTÍCULA
                np=np+1

                !SE COLOCAN LAS 125 PARTICULAS EN LOS VÉRTICES
                rx(np)=0.d00+dble(i)*dble(pa)
                ry(np)=0.d00+dble(j)*dble(pa)
                rz(np)=0.d00+dble(k)*dble(pa)
            
                !SE COLOCA 1 PARTÍCULA POR CADA UNA DE LAS CARAS 3 CARAS
                
                !CARA XY
                rx(np+k3) = rx(np) + pma
                ry(np+k3) = ry(np) + pma
                rz(np+k3) = rz(np)
                !CARA XZ
                rx(np+2*k3) = rx(np) + pma       
                ry(np+2*k3) = ry(np)     
                rz(np+2*k3) = rz(np) + pma
                !CARA YZ
                rx(np+3*k3) = rx(np)
                ry(np+3*k3) = ry(np) + pma            
                rz(np+3*k3) = rz(np) + pma
            end do
        end do
    end do

!   INICIALIZAMOS EL GENERADOR DE NUMEROS ALEATORIOS

    write(*,*)'Escribe un numero entero'
    read(*,*) iflag
    px=mi_random(iflag)

!   INTRODUCIMOS UNA PEQUEÑA MODIFICACION A LAS POSICIONES MEDIANTE NUMEROS ALEATORIOS

    do l=1,npmax
        rx(l) = rx(l) + (2.d00*mi_random(iflag) - 1.d00) * pl/100.d00
        ry(l) = ry(l) + (2.d00*mi_random(iflag) - 1.d00) * pl/100.d00
        rz(l) = rz(l) + (2.d00*mi_random(iflag) - 1.d00) * pl/100.d00
    end do
    

!   LLAMAMOS A LA SUBRUTINA QUE CALCULA LA ENERGÍA POTENCIAL

    call potlj (np,rx,ry,rz,ax,ay,az,epot,dfiv,d2fiv)


!   INTRODUCIMOS VALORES ALEATORIOS DE LA VELOCIDAD (ENTRE -1 Y +1)

    do l = 1,npmax
        vx(l) = 2.d00 * mi_random(iflag) - 1.d00
        vy(l) = 2.d00 * mi_random(iflag) - 1.d00
        vz(l) = 2.d00 * mi_random(iflag) - 1.d00
    end do

!   CALCULAMOS LA ENERGÍA CINETICA QUE NOS DEVUELVE

    px = sum(vx)
    py = sum(vy)
    pz = sum(vz)

    ecin = 1.d00/2.d00 * (sum(vx*vx) + sum(vy*vy) + sum(vz*vz))

! FIJAMOS UN VALOR DE LA ENERGIA TOTAL QUE QUEREMOS, Y OBTENEMOS A PARTIR DE ELLA LA ENERGIA CINETICA

    write(*,*) 'momentos=', px,py,pz
    write(*,*) 'ecin=', ecin,' epot=', epot,' etot=', ecin+epot
    write(*,*) 'que energia total quiere fijar?'
    read(*,*) etotq
    ecinq = etotq - epot

    if (ecinq.le.0.d00) then
        write(*,*) 'error, energia cinetica negativa'
        stop
    endif


!   OBTENEMOS DE NUEVO LAS VELOCIDADES

    factor = dsqrt(ecinq/ecin)
    vx = factor*vx
    vy = factor*vy
    vz = factor*vz

!   RECALCULAMOS PARA VERIFICAR

    px = sum(vx)
    py = sum(vy)
    pz = sum(vz)

    ecin = 0.5d00 * (sum(vx*vx) + sum(vy*vy) + sum(vz*vz))

    write(*,*) 'momentos=', px,py,pz
    write(*,*) 'ecin=', ecin,' epot=', epot,' etot=', ecin+epot



!   PROCEDO A GRABAR FICHEROS Y TERMINO
        
    write(*,*)'Escribe el nombre del archivo de parametros'
    read(*,*) fname 
    write(*,*)'Escribe el nombre del archivo de posiciones, velocidades y aceleraciones: '
    read(*,*) gname

    fname = trim(fname)
    gname = trim(gname)

    open (10,file=fname)
        write(10,*) 'np= ',np*(numk-1),' pl = ',pl,' pli = ',pli,' rc = ',rc,' rc2 = ',rc2
        write(10,*) 'vol = ',vol,' dens = ',dens
        write(10,*) 'etot = ',ecin+epot,' ecin = ',ecin,' epot = ',epot
        write(10,*) 'dfiv = ',dfiv,' d2fiv = ',d2fiv
        write(10,9000) 'fname = ',fname
        write(10,9000) 'gname = ',gname
    close(10)

9000 FORMAT (a24)
    
    open(20,file=gname,form='unformatted')
        write(20) rx,ry,rz,vx,vy,vz,ax,ay,az
    close(20)

end program