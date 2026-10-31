c ======================================================================
c ======================================================================
c === programa modelo2a
c === resolve equacao de crescimento exponencial (dU/dt=gU)
c === progressivo, RK2 e RK4 e Leap-frog
c ===      
c === desenvolvido para grupo de estudos de modelagem do LEPOC      
c === desenvolvido por Afonso Paiva, marco de 2009
c ===
c ======================================================================
c ======================================================================
c === lista de variaveis
c        u1,u2 - numero de individuos (esquema progressivo)
c        v1,v2 - numero de individuos (esquema Runge-Kutta ordem 2)
c        w1,w2 - numero de individuos (esquema Runge-Kutta ordem 4)
c     y0,y1,y2 - numero de individuos (esquema leap-frog)
c     a     - solucao analitica
c     g     - taxa de crescimento
c     dt    - passo de tempo
c     diag  - intervalo de tempo de integracao
c     ttot  - tempo total de integracao
c     tempo - tempo a cada passo de integracao      
c     n     - passo de tempo
c     nn    - numero total de passos de tempo
c ======================================================================
c ======================================================================   
c ===
c === obs: bom grupo de valores para comparar prog com RK2 e RK4
c === g=0.05, ttot=300, diag=10, variar dt = 10, 2, 0.1
c ===
c === obs: bom grupo de valores para comparar leap-frog com os demais
c ===      g=0.05, ttot=300, diag=1, dt = 1
c ===      rodar com u2=u1 e sem filtro de asselin
c ===      plotar curvas entre 0<t>30 e 2<u>6 para ver escada
c ===      Repetir usando filtro
c ===      Repetir u2 a partir de esquema progressivo
c ======================================================================   
c ======================================================================   
      program modelo2

c --- declara variaveis
      real g,tempo,dt,ttot,diag,a,a0,nu
      real u1,u2,v1,v2,w1,w2,y0,y1,y2
      real kv1,kv2,kw1,kw2,kw3,kw4
      integer n,nn

c --- define parametros de integracao numerica
c     obs: ttot e diag devem ser multiplo inteiro de dt      
      ttot = 300.    ! define tempo total de integracao
c      diag = 1.     ! define intervalo para diagnostico dos resultados
c      dt = 1.       ! define passo de tempo
      
      diag = 10.     ! define intervalo para diagnostico dos resultados
      dt = 10.       ! define passo de tempo      
      
      
      nn = anint(ttot/dt)
      write(*,*) ' tempo total de integracao = ',ttot
      write(*,*) ' passo de tempo de integracao = ',dt
      write(*,*) ' numero total de passos de tempo = ',nn
      write(*,*) ' intervalo de tempo para diagnose = ',diag

c --- define parametros do modelo
      g = 0.05 ! taxa de crescimento

c --- define coeficiente do filtro de asselin
      nu = 0.2
      
c --- define condicoes iniciais      
      ci = 2.   ! define condicao inicial
      u1 = ci   ! condicao inicial - esquema progressivo
      v1 = ci   ! condicao inicial - esquema Runge-Kutta ordem 2
      w1 = ci   ! condicao inicial - esquema Runge-Kutta ordem 4
      y1 = ci   ! condicao inicial - esquema leap-frog
      a0 = ci   ! condicao inicial analitica 

c --- abre arquivo de dados e salva condicoes iniciais
      open(1,file='modelo2a.dat',status='unknown')

c --- inicia loop temporal
c     obs: no passo [n], calculo variavel em [n+1] mas salvo em [n]
      do n = 1,nn+1

c ..... integra a equacao no tempo

        u2 = u1 + dt*g*u1 ! (progressivo no tempo)

        kv1 = dt * v1*g       ! (Runge-Kutta ordem 2)
        kv2 = dt * (v1+kv1/2) * g
        v2 = v1 + kv2

        kw1 = dt * ( w1       ) * g ! (Runge-Kutta ordem 4)
        kw2 = dt * ( w1+kw1/2 ) * g
        kw3 = dt * ( w1+kw2/2 ) * g
        kw4 = dt * ( w1+kw3   ) * g
        w2  = w1 + (kw1+2*kw2+2*kw3+kw4)/6

        if (n.eq.1) then   ! leap-frog
          y2 = y1           ! passo inicial - repete tempo inicial
c         y2 = y1 + dt*g*y1 ! passo inicial - progressivo no tempo
        else
          y2 = y0 + 2*dt*g*y1 
          y1 = (1-nu)*y1 + nu*(y0+y2)/2 ! filtro de asselin
        endif
            
c ..... determina o tempo referente a este passo de integracao
c       n=1 indica tempo=0
        tempo = (n-1)*dt

c ..... calcula solucao analitica
        a = a0*exp(g*tempo)

c ..... salva diagnostico
        aux=amod(tempo,diag)
        
c        FMT = "(I4,6(F14.6))"
c        WRITE(*,FMT) ......

        if(aux.eq.0.) write(1,102)n,tempo,u1,v1,w1,y1,a
102     format(I4,6(F18.6))

c ..... resseta variaveis
        u1 = u2
        v1 = v2
        w1 = w2
        y0 = y1
        y1 = y2

c     finaliza loop temporal
      enddo

      close(1)
      stop 'fim normal'
      end
