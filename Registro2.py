#####7segments
import RPi.GPIO as GPIO
import dht11
import time
GPIO.setmode(GPIO.BCM)

## 23 por 16

TEMP = 19
BOTON = 26

GPIO.setup(TEMP, GPIO.IN)
GPIO.setup(BOTON, GPIO.IN)

## 23 por 16
segments = (11,4,16,8,7,10,18,25)
for segment in segments:
    GPIO.setup(segment, GPIO.OUT)
    GPIO.output(segment,0)

## Cambiando 24 por 12 funciona por alguna razon.
## Creo que esta teniendo interferencia con alguna otra cosa
digits=(22,27,17,12)

for digit in digits:
    GPIO.setup(digit,GPIO.OUT)
    GPIO.output(digit,1)

num = {' ':(0,0,0,0,0,0,0),
       '0':(1,1,1,1,1,1,0),
       '1':(0,1,1,0,0,0,0),
       '2':(1,1,0,1,1,0,1),
       '3':(1,1,1,1,0,0,1),
       '4':(0,1,1,0,0,1,1),
       '5':(1,0,1,1,0,1,1),
       '6':(1,0,1,1,1,1,1),
       '7':(1,1,1,0,0,0,0),
       '8':(1,1,1,1,1,1,1),
       '9':(1,1,1,1,0,1,1)}

instance = dht11.DHT11(pin=TEMP)

frecuencia = 2
t_0 = time.time()
result = instance.read()
n = str(result.temperature) + str(result.humidity)
s = str(n).rjust(4)

modo = 1

pulsado = 0
boton = 0


f = open("Data.txt", "w+")
try:
    while True:
	t_1 = time.time()
	# tiempo suficiente como para volver a medir temperatura
	if t_1 - t_0 > frecuencia:
            t_0 = t_1
	    result = instance.read()
            if result.is_valid():
		tiempo = time.ctime()
		f.write(tiempo + "\n")
                n = str(result.temperature) +  str(result.humidity)
                s = str(n).rjust(4)

        # Cambiar modo con el boton
        boton = GPIO.input(26)
        if pulsado == 0 and boton == 1:
            time.sleep(0.1)
            modo = (modo + 1) % 2
            pulsado = 1
        elif boton == 0:
	    pulsado = 0

	# Si esta en modo temp/humedad
        if modo == 0:
            for digit in range(4):
                for loop in range(0,7):
                    GPIO.output(segments[loop],num[s[digit]][loop])                
                if digit == 1:
                    GPIO.output(25,1)
                else:
	            GPIO.output(25,0)
                GPIO.output(digits[digit],0)
                time.sleep(0.005)
                GPIO.output(digits[digit],1)

	# Si esta en modo reloj
        elif modo == 1:
            n = time.ctime()[11:13] + time.ctime()[14:16]
            s = str(n).rjust(4)
            for digit in range(4):
                for loop in range(0,7):
                    GPIO.output(segments[loop], num[s[digit]][loop])
                    if (int(time.ctime()[18:19]) % 2 == 0) and (digit == 1):
                        GPIO.output(25, 1)
                    else:
                        GPIO.output(25, 0)
                GPIO.output(digits[digit], 0)
                time.sleep(0.005)
                GPIO.output(digits[digit], 1)

finally:
    GPIO.cleanup()
    f.close()
