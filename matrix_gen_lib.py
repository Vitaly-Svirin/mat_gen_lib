#matrix_gen_lib
#code for processing genetic matrices
#see more about genetic matrices http://arxiv.org/abs/0803.0888v1

class matrix:
    def __init__(self,array):
        if all([type(array[i])==int or type(array[i])==float for i in range(len(array))]):
            for i in range(len(array)):
                array[i]=[array[i]]
            self.__init__(array)
        elif all([len(array[i])==len(array[i+1]) for i in range(len(array)-1)]):
            self.array=array
            rows=len(array)
            columns=len(array[0])
            self.size=[rows,columns]
        
        else:
            print 'input error'
        

    def __mul__(self,other):
        if type(other)==int or type(other)==float:
            buf=list()
            for i in range(self.size[0]):
                subbuf=list()
                for j in range(self.size[1]):
                    
                    num=self.array[i][j]*other
                    subbuf.append(num)

                buf.append(subbuf)

            x=matrix(buf)
            return x
        if self.size[1]==other.size[0]:
            buf=list()
            for i in range(self.size[0]):
                subbuf=list()
                for j in range(other.size[1]):
                    num=0
                    for k in range(self.size[1]):
                        num+=self.array[i][k]*other.array[k][j]

                    subbuf.append(num)

                buf.append(subbuf)

            x=matrix(buf)
            return x
        else:
            print 'size error'

    def T(self):
        array=[[0 for j in range(self.size[0])] for i in range(self.size[1])]
        for i in range(self.size[0]):
            for j in range(self.size[1]):
                array[j][i]=self.array[i][j]

        x=matrix(array)
        return x

    def c(self):
        array0=self.array
        ln=len(array0)
        array1=list()
        for i in range(ln):
            buf=list()
            for j in range(ln):
                buf.append(array0[i^(ln-1)][j^(ln-1)])

            array1.append(buf)

        x=matrix(array1)
        return x

    def __neg__(self):
        x=self*(-1)
        return x


    def __eq__(self,other):

        r = self.array == other.array
        return r
        
    def __add__(self,other):
        if self.size==other.size:
            buf=list()
            for i in range(self.size[0]):
                subbuf=list()
                for j in range(self.size[1]):
                    num=self.array[i][j]+other.array[i][j]
                    subbuf.append(num)

                buf.append(subbuf)

            x=matrix(buf)
            return x
        else:
            print 'size error'


class multiplet:

#------------G A T C-----------------------------------------------------
    subAlf0=(3,2,2,3)
    subAlf1=(0,0,1,1)
    subAlf2=(0,1,0,1)
    subAlf3=(0,1,1,0)
    alf=[subAlf0,subAlf1,subAlf2,subAlf3]
    
    dic1={subAlf1[0]:'G',subAlf1[1]:'A',subAlf1[2]:'T',subAlf1[3]:'C'}
    dic2={subAlf2[0]:'G',subAlf2[1]:'A',subAlf2[2]:'T',subAlf2[3]:'C'}
    dic3={subAlf3[0]:'G',subAlf3[1]:'A',subAlf3[2]:'T',subAlf3[3]:'C'}
    
    
    dicSk1={'G':subAlf1[0],'A':subAlf1[1],'T':subAlf1[2],'C':subAlf1[3]}
    dicSk2={'G':subAlf2[0],'A':subAlf2[1],'T':subAlf2[2],'C':subAlf2[3]}
    dicSk3={'G':subAlf3[0],'A':subAlf3[1],'T':subAlf3[2],'C':subAlf3[3]}
    dicSkF = lambda self,i:{'G':self.alf[i][0],'A':self.alf[i][1],'T':self.alf[i][2],'C':self.alf[i][3]}
    
    dicDSk123= {'G':subAlf1[0]*4+ subAlf2[0]*2+subAlf3[0], 'A':subAlf1[1]*4+ subAlf2[1]*2+subAlf3[1],
               'T':subAlf1[2]*4+ subAlf2[2]*2+subAlf3[2], 'C':subAlf1[3]*4+ subAlf2[3]*2+subAlf3[3]}

    dicBSk123= {'G':[subAlf1[0],subAlf2[0],subAlf3[0]], 'A':[subAlf1[1],subAlf2[1],subAlf3[1]],
               'T':[subAlf1[2],subAlf2[2],subAlf3[2]], 'C':[subAlf1[3],subAlf2[3],subAlf3[3]]}

    dicBSk12= {'G':[subAlf1[0],subAlf2[0]], 'A':[subAlf1[1],subAlf2[1]],
               'T':[subAlf1[2],subAlf2[2]], 'C':[subAlf1[3],subAlf2[3]]}
    dicBSk13= {'G':[subAlf1[0],subAlf2[0]], 'A':[subAlf1[1],subAlf2[1]],
               'T':[subAlf1[2],subAlf3[2]], 'C':[subAlf1[3],subAlf3[3]]}
    dicBSk21= {'G':[subAlf2[0],subAlf1[0]], 'A':[subAlf2[1],subAlf1[1]],
               'T':[subAlf2[2],subAlf1[2]], 'C':[subAlf2[3],subAlf1[3]]}
    dicBSk23= {'G':[subAlf2[0],subAlf3[0]], 'A':[subAlf2[1],subAlf2[1]],
               'T':[subAlf2[2],subAlf3[2]], 'C':[subAlf2[3],subAlf3[3]]}
    dicBSk31= {'G':[subAlf3[0],subAlf1[0]], 'A':[subAlf3[1],subAlf1[1]],
               'T':[subAlf3[2],subAlf1[2]], 'C':[subAlf3[3],subAlf1[3]]}
    dicBSk32= {'G':[subAlf3[0],subAlf2[0]], 'A':[subAlf3[1],subAlf2[1]],
               'T':[subAlf3[2],subAlf2[2]], 'C':[subAlf3[3],subAlf2[3]]}

    dic12={(subAlf1[0]*2+subAlf2[0]): 'G',(subAlf1[1]*2+subAlf2[1]): 'A',
           (subAlf1[2]*2+subAlf2[2]): 'T',(subAlf1[3]*2+subAlf2[3]): 'C' }
    dic13={(subAlf1[0]*2+subAlf3[0]): 'G',(subAlf1[1]*2+subAlf3[1]): 'A',
           (subAlf1[2]*2+subAlf3[2]): 'T',(subAlf1[3]*2+subAlf3[3]): 'C' }
    dic21={(subAlf2[0]*2+subAlf1[0]): 'G',(subAlf2[1]*2+subAlf1[1]): 'A',
           (subAlf2[2]*2+subAlf1[2]): 'T',(subAlf2[3]*2+subAlf1[3]): 'C' }
    dic23={(subAlf2[0]*2+subAlf3[0]): 'G',(subAlf2[1]*2+subAlf3[1]): 'A',
           (subAlf2[2]*2+subAlf3[2]): 'T',(subAlf2[3]*2+subAlf3[3]): 'C' }
    dic31={(subAlf3[0]*2+subAlf1[0]): 'G',(subAlf3[1]*2+subAlf1[1]): 'A',
           (subAlf3[2]*2+subAlf1[2]): 'T',(subAlf3[3]*2+subAlf1[3]): 'C' }
    dic32={(subAlf3[0]*2+subAlf2[0]): 'G',(subAlf3[1]*2+subAlf2[1]): 'A',
           (subAlf3[2]*2+subAlf2[2]): 'T',(subAlf3[3]*2+subAlf2[3]): 'C' }

    

    
    
    def __init__(self,length,num,AlfA,AlfB):
        
        AlfOrder = AlfA*10 + AlfB        
        self.len=length

        if AlfOrder == 12:
            
            self.dic=multiplet.dic12
            
        if AlfOrder == 13:
            
            self.dic=multiplet.dic13
            
        if AlfOrder == 21:
            
            self.dic=multiplet.dic21

        if AlfOrder == 23:
            
            self.dic=multiplet.dic23

        if AlfOrder == 31:
            
            self.dic=multiplet.dic31

        if AlfOrder == 32:
            
            self.dic=multiplet.dic32

        self.name=''
        num123=[]
        num12=[]
        num13=[]
        num21=[]
        num23=[]
        num31=[]
        num32=[]
        
        
        for i in range(1,self.len+1):
            numL =(num / (2**self.len) ) / (2**(self.len-i))-((num / (2**self.len) ) / (2**(self.len-i+1)))*2  
            numR =(num %(2**self.len)) / (2**(self.len-i))-((num %(2**self.len)) / (2**(self.len-i+1)))*2
            letter=self.dic[numL*2+numR]
            num123.append(multiplet.dicBSk123[letter])
            num12.append(multiplet.dicBSk12[letter])
            num13.append(multiplet.dicBSk13[letter])
            num21.append(multiplet.dicBSk21[letter])
            num23.append(multiplet.dicBSk23[letter])
            num31.append(multiplet.dicBSk31[letter])
            num32.append(multiplet.dicBSk32[letter])
            self.name += letter


        self.num123=0
        self.num12=0
        self.num13=0
        self.num21=0
        self.num23=0
        self.num31=0
        self.num32=0
        
        
        for i in range(self.len):
            self.num123 += num123[i][0]*2**(3*self.len-(i+1)) + num123[i][1]*2**(2*self.len-(i+1)) + num123[i][2]*2**(self.len-(i+1))
            self.num12 += num12[i][0]*2**(2*self.len-(i+1)) + num12[i][1]*2**(self.len-(i+1))
            self.num13 += num13[i][0]*2**(2*self.len-(i+1)) + num13[i][1]*2**(self.len-(i+1))
            self.num21 += num21[i][0]*2**(2*self.len-(i+1)) + num21[i][1]*2**(self.len-(i+1))
            self.num23 += num23[i][0]*2**(2*self.len-(i+1)) + num23[i][1]*2**(self.len-(i+1))
            self.num31 += num31[i][0]*2**(2*self.len-(i+1)) + num31[i][1]*2**(self.len-(i+1))
            self.num32 += num32[i][0]*2**(2*self.len-(i+1)) + num32[i][1]*2**(self.len-(i+1))




    def initByName(self,name):
        self.len = len(name)
        self.name = name
        
        num123=[]
        num12=[]
        num13=[]
        num21=[]
        num23=[]
        num31=[]
        num32=[]
        

        for i in range(self.len):
            num123.append(multiplet.dicBSk123[self.name[i]])
            num12.append(multiplet.dicBSk12[self.name[i]])
            num13.append(multiplet.dicBSk13[self.name[i]])
            num21.append(multiplet.dicBSk21[self.name[i]])
            num23.append(multiplet.dicBSk23[self.name[i]])
            num31.append(multiplet.dicBSk31[self.name[i]])
            num32.append(multiplet.dicBSk32[self.name[i]])
            
        
        self.num123=0
        self.num12=0
        self.num13=0
        self.num21=0
        self.num23=0
        self.num31=0
        self.num32=0
        
        
        for i in range(self.len):
            self.num123 += num123[i][0]*2**(3*self.len-(i+1)) + num123[i][1]*2**(2*self.len-(i+1)) + num123[i][2]*2**(self.len-(i+1))
            self.num12 += num12[i][0]*2**(2*self.len-(i+1)) + num12[i][1]*2**(self.len-(i+1))
            self.num13 += num13[i][0]*2**(2*self.len-(i+1)) + num13[i][1]*2**(self.len-(i+1))
            self.num21 += num21[i][0]*2**(2*self.len-(i+1)) + num21[i][1]*2**(self.len-(i+1))
            self.num23 += num23[i][0]*2**(2*self.len-(i+1)) + num23[i][1]*2**(self.len-(i+1))
            self.num31 += num31[i][0]*2**(2*self.len-(i+1)) + num31[i][1]*2**(self.len-(i+1))
            self.num32 += num32[i][0]*2**(2*self.len-(i+1)) + num32[i][1]*2**(self.len-(i+1))




    def  __mul__ (self,other):
        name=self.name + other.name
        x=multiplet(1,0,1,2)
        x.initByName(name)
        return x
             

    def __len__(self):
        return self.len

    def c (self):
        len=self.len
        num12=self.num12 ^ (2**(2*len)-1)
        x=multiplet(len,num12,1,2)
        return x

    def read (self,seqN,seqA):
        name=self.name
        out=0
        ln=len(seqN)
        for i in range(ln):
            out += self.dicSkF(seqA[i])[name[seqN[i]-1]] * 2**(ln-i-1)
        return out

    def readVec (self,NAlf):
        name=self.name
        out=list()
        for i in range(self.len):
            out.append(self.dicSkF(NAlf)[name[i]])
        return out

    def sign (self,seqN,seqA,rule,sv):
        num=self.read(seqN,seqA)
        if num in rule:
            return sv
        else:
            return -sv



class genomat:
    def __init__(self,arry1,arry2,alfA,alfB):
        if len(arry1)<>len(arry2):
            print 'length of arry1 must be equal to length of arry2'
            return 1
        if len((arry1)*len(arry2))==0:
            print 'arrys must have not zero length'
            return 2

        gm=list()
        lena=len(arry1)
        for i in range(lena):
            str=list()
            for j in range(lena):
                lenm=1
                while 4**lenm < lena**2:
                    lenm += 1
                
                num=arry1[i]*2**(lenm)+arry2[j]
                x=multiplet(lenm,num,alfA,alfB)
                mult=x.name
                str.append(mult)

            gm.append(str)

        self.gm = gm
        self.size=lena

                
    def pgm(self):
        for i in range(len(self.gm)):
            print self.gm[i]


    def read (self,seqN,seqA,p):
        ln=self.size
        numm=list()
        for i in range(ln):
            nums=list()
            for j in range(ln):
                name = self.gm[i][j]
                x=multiplet(1,0,1,2)
                x.initByName(name)
                nums.append(x.read(seqN,seqA))

            numm.append(nums)

        if p==0:        
            return numm
        if p==2:
            for i in range(ln):
                print numm[i]
            return numm
        if p==1:
            for i in range(ln):
                print numm[i]

    def readVec (self,NAlf,p):
        ln=self.size
        numm=list()
        for i in range(ln):
            nums=list()
            for j in range(ln):
                name = self.gm[i][j]
                x=multiplet(1,0,1,2)
                x.initByName(name)
                nums.append(x.readVec(NAlf))

            numm.append(nums)

        if p==0:        
            return numm
        if p==2:
            for i in range(ln):
                print numm[i]
            return numm
        if p==1:
            for i in range(ln):
                print numm[i]

    def c (self):
        ln=self.size
        newm=genomat([0],[0],1,2)
        buf=list()
        for i in range(ln):
            bufs=list()
            for j in range(ln):
                name=self.gm[i][j]
                x=multiplet(1,0,1,2)
                x.initByName(name)
                y=x.c()
                bufs.append(y.name)
            buf.append(bufs)
        newm.gm=buf
        return newm


    def sign (self,seqN,seqA,rule,sv,p):
        ln=self.size
        sign=list()
        for i in range(ln):
            ssign=list()
            for j in range(ln):
                name=self.gm[i][j]
                x=multiplet(1,0,1,2)
                x.initByName(name)
                ssign.append(x.sign(seqN,seqA,rule,sv))
            sign.append(ssign)
        if p==0:        
            return sign
        if p==2:
            for i in range(ln):
                print sign[i]
            return sign
        if p==1:
            for i in range(ln):
                print sign[i]


    def decomp(self,seqNS,seqAS,ruleS,sv,seqND,seqAD,nc):
        sign = self.sign(seqNS,seqAS,ruleS,sv,0)
        return self.dcl(sign,seqND,seqAD,nc)


    def dcl(self,slist,seqND,seqAD,nc):
        sign = slist
        mosaic = self.read(seqND,seqAD,0)
        ln=self.size
        ef=list()
        for i in range(ln):
            ei=list()
            for k in range(ln):
                subei=list()
                for n in range(ln):
                    if mosaic[k][n] == i:
                        subei.append(sign[k][n])

                    else:
                        subei.append(0)
                ei.append(subei)

            ei=matrix(ei)
            ef.append(ei)
        return self.check(ef,nc)

    def check(self,basis,nc):
        ln=self.size
        buf=list()
        numbuf=list()
        null=matrix([[0 for i in xrange(ln)] for j in xrange(ln)])
        E=matrix([[int(not i^j) for i in xrange(ln)] for j in xrange(ln)])
        if nc==0:
            for i in range(ln):
                subbuf=list()
                subnumbuf=list()
                for j in range(ln):
                    e=basis[i]*basis[j]
                    if e in basis:
                        subbuf.append('B'+str(basis.index(e)))
                        subnumbuf.append(1)
                    elif -e in basis:
                        subbuf.append('-B'+str(basis.index(-e)))
                        subnumbuf.append(-1)
                    elif e == null:
                        subbuf.append('0')
                        subnumbuf.append(0)
##                      print 'no algebra'
                    else:
                        return 0

            
##              print subbuf            
                buf.append(subbuf)
                numbuf.append(subnumbuf)

            r=[basis,buf,numbuf]
            
        else:
            for n in range(ln):
                subbuf=list()
                subnumbuf=list()
                for k in range(ln):
                    e=matrix([[0 for i in range(ln)] for i in range(ln)])
                    for i in range(ln):
                        e.array[i][i^n^k]=basis[n].array[i][i^n]*basis[k].array[i^n][i^n^k]
                    if e in basis:
                        if e ==E:
                            subbuf.append('1')
                        elif e==null:
                            subbuf.append('0')
                        else:
                            subbuf.append('B'+str(basis.index(e)))
                        subnumbuf.append(1)
                    elif -e in basis:
                        if e ==E:
                            subbuf.append('-1')
                        else:
                            subbuf.append('-B'+str(basis.index(-e)))
                        subnumbuf.append(-1)
                    elif e==null:
                        subbuf.append('0')
                        subnumbuf.append(0)
                    else:
##                              print 'no algebra'
                        return 0

            
##                      print subbuf            
                buf.append(subbuf)
                numbuf.append(subnumbuf)

            r=[basis,buf,numbuf]
        return r