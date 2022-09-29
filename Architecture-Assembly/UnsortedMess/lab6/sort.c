/**
 * Sort an array using O(n^2) bubble sort.
 * @param data  array of 64-bit integers
 * @param count  number of elements in the input array
 */
void sort(long *data, long count)
{
   for (long last = count - 1; last > 0; --last) {
       for (long i = 0; i < last; ++i) {
           if (data[i + 1] < data[i]) {
               long t = data[i + 1];
               data[i + 1] = data[i];
               data[i] = t;
           }
       }
   }
}
