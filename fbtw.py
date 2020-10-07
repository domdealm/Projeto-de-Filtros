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
        if self.tipo == 'Passa-baixa':
            n = int(np.ceil(
                np.log10((np.power(10,(abs(self.Ap)/10))-1)/(np.power(10,(abs(self.As)/10))-1))/(2*np.log10(self.Wp/self.Ws))
            ))
        if self.tipo == "Passa-alta":
            n = int(np.ceil(
                np.log10((np.power(10,(abs(self.Ap)/10))-1)/(np.power(10,(abs(self.As)/10))-1))/(2*np.log10(self.Ws/self.Wp))
            ))
        if self.tipo == "Passa-faixa":
            #aproximando a ordem pelo pior caso
            n1 = int(np.ceil(
                np.log10((np.power(10,(abs(self.Ap)/10))-1)/(np.power(10,(abs(self.As)/10))-1))/(2*np.log10(self.Wp2/self.Ws2))
            ))
            n2 = int(np.ceil(
                np.log10((np.power(10,(abs(self.Ap)/10))-1)/(np.power(10,(abs(self.As)/10))-1))/(2*np.log10(self.Ws1/self.Wp1))
            ))
            if(n1>n2):
                n = n1
            else:
                n = n2
        if self.tipo == "Rejeita-faixa":
            #aproximando a ordem pelo pior caso
            n1 = int(np.ceil(
                np.log10((np.power(10,(abs(self.Ap)/10))-1)/(np.power(10,(abs(self.As)/10))-1))/(2*np.log10(self.Wp1/self.Ws1))
            ))
            n2 = int(np.ceil(
                np.log10((np.power(10,(abs(self.Ap)/10))-1)/(np.power(10,(abs(self.As)/10))-1))/(2*np.log10(self.Ws2/self.Wp2))
            ))
            if(n1>n2):
                n = n1
            else:
                n = n2
        return n
    
    def freq_corte(self):
        wc = None
        wc1 = 0
        wc2 = 0
        n = self.ordem()
        if self.tipo == 'Passa-baixa':
            wc = self.Wp/(np.power((np.power(10,(abs(self.Ap)/10))-1),(1/(2*n))))
        if self.tipo == 'Passa-alta':
            wc = self.Wp*(np.power((np.power(10,(abs(self.Ap)/10))-1),(1/(2*n))))
        if self.tipo == 'Passa-faixa':
            wc1 = self.Wp1*(np.power((np.power(10,(abs(self.Ap)/10))-1),(1/(2*n))))
            wc2 = self.Wp2/(np.power((np.power(10,(abs(self.Ap)/10))-1),(1/(2*n))))
            wc = np.sqrt(wc1*wc2)
        if self.tipo == 'Rejeita-faixa':
            wc1 = self.Wp1/(np.power((np.power(10,(abs(self.Ap)/10))-1),(1/(2*n))))
            wc2 = self.Wp2*(np.power((np.power(10,(abs(self.Ap)/10))-1),(1/(2*n))))
            wc = np.sqrt(wc1*wc2)
        wout = [wc,wc1,wc2]
        return wout
    
    def raizes(self):
        n = self.ordem()
        sk = np.zeros(n, dtype=complex)
        for i in range(1,n+1):
            sk[i-1] = -np.sin((np.pi*((2*i)-1))/(2*n))+(1j*np.cos((np.pi*((2*i)-1))/(2*n)))
        return sk
    
    def fun_trans(self):
        k = self.raizes()
        wc = self.freq_corte()
        Bw = (wc[2]-wc[1])
        self.nden = np.real(np.poly(k)) #denominador normalizado
        den = np.zeros(len(k)) #criando o array para o denominador transformado
        if self.tipo == 'Passa-baixa':
            [num,den] = sig.lp2lp(self.nden[-1],self.nden,wc[0])
            FT = sig.TransferFunction(num,den) #gerando a FT para a freq. de corte
        if self.tipo == 'Passa-alta':
            [num,den] = sig.lp2hp(self.nden[-1],self.nden,wc[0])
            FT=sig.TransferFunction(num,den)
        if self.tipo == 'Passa-faixa':
            [num,den] = sig.lp2bp(self.nden[-1],self.nden,wc[0],Bw)
            FT=sig.TransferFunction(num,den)
        if self.tipo == 'Rejeita-faixa':
            [num,den] = sig.lp2bs(self.nden[-1],self.nden,wc[0],Bw)
            FT=sig.TransferFunction(num,den)
        return FT
    


    




