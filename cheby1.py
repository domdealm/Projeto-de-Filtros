import numpy as np
from scipy import signal as sig
#Tipos de filtro:
#PB = Passa baixas
#PA = Passa altas
#PF = Passa faixa
#RF = Rejeita faixa
class filtro_cb1:
   
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
        if self.tipo == 'Passa-baixa':
            n = int(np.ceil(np.arccosh((np.sqrt(np.power(10,abs(self.As)/10)-1))/(np.sqrt(np.power(10,abs(self.Ap)/10)-1)))/(np.arccosh(self.Ws/self.Wp))))
        if self.tipo == "Passa-alta":
            n = int(np.ceil(np.arccosh((np.sqrt(np.power(10,abs(self.As)/10)-1))/(np.sqrt(np.power(10,abs(self.Ap)/10)-1)))/(np.arccosh(self.Wp/self.Ws))))
        if self.tipo == "Passa-faixa":
            n = int(np.ceil(np.arccosh((np.sqrt(np.power(10,abs(self.As)/10)-1))/(np.sqrt(np.power(10,abs(self.Ap)/10)-1)))/(np.arccosh((self.Ws2 - self.Ws1)/(self.Wp2 - self.Wp1)))))
        if self.tipo == "Rejeita-faixa":
            n = int(np.ceil(np.arccosh((np.sqrt(np.power(10,abs(self.As)/10)-1))/(np.sqrt(np.power(10,abs(self.Ap)/10)-1)))/(np.arccosh((self.Wp2 - self.Wp1)/(self.Ws2 - self.Ws1)))))
        return n
    
    def epso(self):
        epsi = np.sqrt(np.power(10,(abs(self.Ap)/10)) - 1)
        return epsi

    def freq_corte(self):
        wc = None
        n = self.ordem()
        eps = self.epso()
        if self.tipo == 'Passa-baixa':
            wc = self.Wp*np.cosh((1/n)*np.arccosh(1/eps))
        if self.tipo == 'Passa-alta':
            wc = self.Wp/np.cosh((1/n)*np.arccosh(1/eps))
        if self.tipo == 'Passa-faixa':
            #aproximando por Wp
            wc = np.sqrt(self.Wp1*self.Wp2)
        if self.tipo == 'Rejeita-faixa':
            #aproximando por Wp
            wc = np.sqrt(self.Wp1*self.Wp2)
        return wc
    
    def raizes(self):
        n = self.ordem()
        eps = self.epso()
        ak = np.zeros(n, dtype=complex)
        wk = np.zeros(n, dtype=complex)
        for i in range(1, n + 1):
            ak[i-1] = -1*np.sinh((1/n)*np.arcsinh(1/eps))*np.sin( (np.pi/(2*n))*(2*i - 1))
            wk[i-1] = 1*np.cosh((1/n)*np.arcsinh(1/eps))*np.cos((np.pi/(2*n))*(2*i - 1))
        raizes = ak + 1j*wk
        return raizes
        
    
    def fun_trans(self):
        raizes = self.raizes()
        wc = self.freq_corte()
        n = self.ordem()
        eps = self.epso()
        nden = np.real(np.poly(raizes))
        if(n%2==0):
            nnum = nden[-1]*(1/np.sqrt(1+np.power(eps,2)))
        else:
            nnum = nden[-1]
        if self.tipo == 'Passa-baixa':
            [num,den] = sig.lp2lp(nnum,nden,wc)
        if self.tipo == 'Passa-alta':
            [num,den] = sig.lp2hp(nnum,nden,wc)
        if self.tipo == 'Passa-faixa':
            Bw = self.Wp2 - self.Wp1
            [num,den] = sig.lp2bp(nnum,nden,wc,Bw)
        if self.tipo == 'Rejeita-faixa':
            Bw = self.Wp2 - self.Wp1
            [num,den] = sig.lp2bs(nnum,nden,wc,Bw)
        FT = sig.TransferFunction(num,den)
        return FT
    


    




