unsigned long factorial(unsigned long n){
  if(n == 0) return 1;
  return n*factorial(n-1);
}

unsigned long factorialWhile(unsigned long n){
  long result = 1;
  while(n > 1){
    result *= n;
    n = n - 1;
  }
  return result;
}


long inc_a(long a) {return a++;}
long inc_b(long b) {return ++b;}
