import numpy as np
from scipy import signal as sig
#Tipos de filtro:
#PB = Passa baixas
#PA = Passa altas
#PF = Passa faixa
#RF = Rejeita faixa
class filtro_bw:
   
    def __init__(self,tipo,**kwargs):
        self.tipo = tipo
        self.Ap = kwargs['Ap']          
        self.As = kwargs['As']          
        #Parâmetros para Passa-baixa/Passa-alta
        self.Ws = kwargs.get('Ws',0) #quando não tiver valor atribui 0 por padrão
        self.Wp = kwargs.get('Wp',0)
        #Parâmetros para Rejeita-faixa/Passa-faixa
        self.Ws1 = kwargs.get('Ws1',0)  
        self.Ws2 = kwargs.get('Ws2',0)
        self.Wp1 = kwargs.get('Wp1',0)
        self.Wp2 = kwargs.get('Wp2',0)

    def ordem(self):
        n = 0
        if self.tipo == 'PB':
            n = int(np.ceil(np.log10((np.power(10,((-self.Ap/10)-1))-1)/(np.power(10,((-self.As/10)-1))-1))/2*np.log10(self.Wp/self.Ws)))
        if self.tipo == "PA":
            n = int(np.ceil(np.log10((np.power(10,((-self.Ap/10)-1))-1)/(np.power(10,((-self.As/10)-1))-1))/2*np.log10(self.Ws/self.Wp)))
        return n
    
    def freq_corte(self):
        wc = None
        n = self.ordem()
        if self.tipo == 'PB':
            wc = self.Wp/(np.power(np.power(10,((-self.Ap/10)-1))-1,1/(2*n)))
        if self.tipo == 'PA':
            wc = self.Wp*(np.power(np.power(10,((-self.Ap/10)-1))-1,1/(2*n)))
        return wc
    
    def raizes(self):
        n = self.ordem()
        sk = np.zeros(n, dtype=complex)
        for i in range(1,n+1):
            sk[i-1] = -np.sin((np.pi*((2*i)-1))/(2*n))+(1j*np.cos((np.pi*((2*i)-1))/(2*n)))
        return sk
    
    def fun_trans(self):
        k = self.raizes()
        n = self.ordem()
        wc = self.freq_corte()
        self.nden = np.real(np.poly(k)) #denominador normalizado
        den = np.zeros(n) #criando o array para o denominador transformado
        for i in range(0,n):
            den[i] = self.nden[i] * np.power(wc,i)
        num = den[-1]
        FT = sig.TransferFunction(num,den) #gerando a FT para a freq. de corte
        return FT


    




