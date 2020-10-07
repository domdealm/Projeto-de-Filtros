from fbtw import *
from cheby1 import *
import inquirer as iq
from matplotlib import pyplot as plt
from scipy import signal 
print("#######################\n\n")
print("Cálculo de filtros\n")
print("\n\n#######################")
############ menu 1 ##############
men1 = [iq.List('M1',message="Qual filtro você deseja calcular?",
        choices=[('Butterworth',1),('Chebyshev 1',2),('Chebyshev 2',3)],carousel=True),]
opt1 = iq.prompt(men1)
############ menu 2 ##############
men2 = [iq.List('M2',message="Escolha o tipo do filtro",
        choices=[('Passa-baixa','Passa-baixa'),('Passa-alta','Passa-alta'),('Passa-faixa','Passa-faixa'),('Rejeita-faixa','Rejeita-faixa')],carousel=True),]
opt2 = iq.prompt(men2) 
############ Butterworth selecionado ##############
if(opt1.get('M1')==1):

    if(opt2.get('M2')=='Passa-baixa' or opt2.get=='Passa-alta'):
        apdb = int(input("Digite o valor de Ap:"))
        asdb = int(input("Digite o valor de As:"))
        wp = int(input("Digite o valor de Wp:"))
        ws = int(input("Digite o valor de Ws:"))
        tipo = opt2.get('M2')
        filtro = filtro_bw(tipo,Ap = apdb,As = asdb,Ws = ws,Wp = wp)
    else:
        apdb = int(input("Digite o valor de Ap:"))
        asdb = int(input("Digite o valor de As:"))
        wp1 = int(input("Digite o valor de Wp1:"))
        wp2 = int(input("Digite o valor de Wp2:"))
        ws1 = int(input("Digite o valor de Ws1:"))
        ws2 = int(input("Digite o valor de Ws2:"))
        tipo = opt2.get('M2')
        filtro = filtro_bw(tipo,Ap = apdb,As = asdb,Ws1 = ws1,Wp1 = wp1,Ws2 = ws2,Wp2 = wp2)

    fft = filtro.fun_trans()
    w,mag,phase = signal.bode(fft,n=400)
    plt.figure()
    plt.semilogx(w, mag)    # Bode magnitude plot
    plt.title('Filtro Butterworth ('+str(opt2.get('M2'))+')')
    plt.xlabel('Frequência')
    plt.ylabel('Ganho (dB)')
    #plt.xlim(0,1e5)
    plt.ylim(-50,1)
    plt.grid(True)
    plt.show()
        
############ Chebyshev 1 selecionado ##############
        
if(opt1.get('M1')==2):

    if(opt2.get('M2')=='PB' or opt2.get=='PA'):
        apdb = int(input("Digite o valor de Ap:"))
        asdb = int(input("Digite o valor de As:"))
        wp = int(input("Digite o valor de Wp:"))
        ws = int(input("Digite o valor de Ws:"))
        tipo = opt2.get('M2')
        filtro = filtro_cb1(tipo,Ap = apdb,As = asdb,Ws = ws,Wp = wp)
    else:
        apdb = int(input("Digite o valor de Ap:"))
        asdb = int(input("Digite o valor de As:"))
        wp1 = int(input("Digite o valor de Wp1:"))
        wp2 = int(input("Digite o valor de Wp2:"))
        ws1 = int(input("Digite o valor de Ws1:"))
        ws2 = int(input("Digite o valor de Ws2:"))
        tipo = opt2.get('M2')
        filtro = filtro_cb1(tipo,Ap = apdb,As = asdb,Ws1 = ws1,Wp1 = wp1,Ws2 = ws2,Wp2 = wp2)
    
    fft = filtro.fun_trans()
    w,mag,phase = signal.bode(fft,n=400)
    plt.figure()
    plt.semilogx(w, mag)    # Bode magnitude plot
    plt.title('Filtro Chebyshev 1 ('+str(opt2.get('M2'))+')')
    plt.xlabel('Frequência')
    plt.ylabel('Ganho (dB)')
    #plt.xlim(0,1e5)
    plt.ylim(-50,1)
    plt.grid(True)
    plt.show()

else:
    print("Em desenvolvimento 2")
