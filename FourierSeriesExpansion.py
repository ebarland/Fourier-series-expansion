from numpy import *
import matplotlib.pyplot as plt

############ FOURIER SERIES EXPANSION ############
# BASES ITSELF ON THE SAME PARAMETERS AS SQUARE WAVE GENERATOR
# TO MAKE A FOURIER SERIES EXPANSION THAT MATCHES THE SQUARE WAVE
def fourier_series_expansion(t1, t2, C1, C2, D, R, N, t):
   
    T = float(t2-t1)
    L = float(T/2) #The bounds of integration will be -L to L
    
    a0 = 1.0/(2.0*L)*(C1*L+C2*L) + D #Calculate a0 and add the vertical shift 
    
    an = zeros(N)
    for n in range(1,N):
        an[n] = 1.0/(n*pi)*(C2*sin(n*pi)-C1*sin(n*pi))
    
    bn = zeros(N)
    for n in range(1,N): 
        bn[n] = 1.0/(n*pi)*(C2*(cos(-n*pi)-1)+C1*(1-cos(n*pi)))
        
    fourier_series_wave = zeros(R)
    
    for i in range(0, R):
        fourier_series_wave[i] += a0 #Add the constant value a0 to the current time step 
        for n in range(0, N): 
            fourier_series_wave[i] += an[n]*cos(pi*n*t[i]/L) + bn[n]*sin(pi*n*t[i]/L) #Adding up all the terms of the fourier transform for the current time step
    
    return fourier_series_wave
    
############ SQUARE WAVE GENERATOR ############
# USES SIN WAVE TO ESTIMATE IF TOP OR BOTTOM LINE WILL BE DRAWN
# SIN VALUE ABOVE VERTICAL MID POINT ADDS HIGH CONSTANT VALUE
# SIN VALUE BELOW VERTICAL MID POINT ADDS LOW CONSTANT VALUE
def square_wave_generator(t1, t2, C1, C2, D, R, t):    
    A = (C1-C2)/2.0 #Calculate the amplitude of the wave
    B = 2.0*pi/(t2-t1) #B is equal to 2Pi/period
    square_wave = zeros(R) #Generate empty matrix of resolution length
    
    for i in range(0, R):
        if A*sin(B*t[i]) + D >= D:
            square_wave[i] = C1 + D
        else:
            square_wave[i] = C2 + D
    return square_wave    

############ FILTERS ############
#Only lets low frequencies pass
def low_pass_filter(signal, dt, tau):
    a = dt/(dt+tau)
    filtered_signal = zeros(len(signal))
    filtered_signal=signal;
    for i in range(1, len(signal)):
       filtered_signal[i] = a*filtered_signal[i] + (1-a)*filtered_signal[i-1]
    return filtered_signal
    
#Only lets high frequency signal through. 
def high_pass_filter(signal, dt, tau):
    a = tau/(dt + tau)
    filtered_signal = zeros(len(signal))
    filtered_signal[0] = signal[0];
    for i in range(1, len(signal)):
       filtered_signal[i] = a*filtered_signal[i-1] + a*(signal[i] - signal[i-1])
    return filtered_signal
    
############ PARAMETERS ############

# START AND END OF PERIOD
ta = 0 #Start of signal pulse
tb = pi #End of signal pulse

# VERTICAL SHIFT 
vshift = 0

# HIGH AND LOW CONSTANT VALUES OF SQUARE WAVE
V1 = 1 #The highest value of the square wave
V2 = -1 #The lowest value of the square wave

# LENGTH OF TIMELINE
tStart = 0.0 #Start of the plot timeline
tEnd = 10.0 #End of the plot timeline

# RESOLUTION OF TIMELINE AND FOURIER TRANSFORM SIGMA
Rt = 10000 #Resolution of time-line
Rf = 20 #Resolution of fourier series

dt = (tEnd - tStart)/Rt #Timestep

# GENERATED TIMELINE
time = linspace(tStart, tEnd, Rt, endpoint = False)

############ PLOTTING AND FUNCTION CALLS ############
# CALL FUNCTIONS TO GENERATE FOURIER SERIES WAVE AND SQUARE WAVE
fourier_series_wave = fourier_series_expansion(ta, tb, V1, V2, vshift, Rt, Rf, time)
square_wave = square_wave_generator(ta, tb, V1, V2, vshift, Rt, time)

#PLOT DATA SETS 
plt.subplot(311)
plt.plot(time, fourier_series_wave, label = 'Fourier expansion')
plt.plot(time, square_wave, label = 'Real square wave')
plt.xlabel('Time [t]')
plt.ylabel('Amplitude [A(t)]')
plt.tight_layout()
plt.legend()

plt.subplot(312)
plt.plot(time, high_pass_filter(fourier_series_wave, dt, 0.07), label = 'High pass filtered')
plt.xlabel('Time [t]')
plt.ylabel('Amplitude [A(t)]')
plt.legend()

plt.subplot(313)
plt.plot(time, low_pass_filter(fourier_series_wave, dt, 0.07), label = 'Low pass filtered')
plt.xlabel('Time [t]')
plt.ylabel('Amplitude [A(t)]')
plt.legend()

