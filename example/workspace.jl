using ChaoticCommunications

generator = Generator(100) 
modulator = CSKModulator([Logistic(), Logistic()], 100)
modulator(generator.bits)

