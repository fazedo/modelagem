import time
import calor
import calor2
import calor3

start_time = time.time()
[x, u] = calor2.calor(1.0, 0.1, 201, 2001)
elapsed_time = time.time() - start_time
print('elapsed_time:', elapsed_time)
#print(u)

start_time = time.time()
calor.calor(1.0, 0.1, 201, 2001)
elapsed_time = time.time() - start_time

print('elapsed_time:', elapsed_time)
#print(x)
#print(u)

start_time = time.time()
[x, u] = calor3.calor(1.0, 0.1, 201, 2001)
elapsed_time = time.time() - start_time
print('elapsed_time:', elapsed_time)
print(u)
