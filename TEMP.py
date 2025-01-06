import serial
import time

# Replace 'COM1' with the correct port name, e.g., '/dev/ttyUSB0' on Linux
ser = serial.Serial('COM3', 9600, timeout=1)

print("Listening on serial port...")
time.sleep(2)
command = "C2\n"
ser.write(command.encode('utf-8'))
while True:

    data = ser.read(1)  # Read one byte at a time
    if data:
        print(data.hex())
