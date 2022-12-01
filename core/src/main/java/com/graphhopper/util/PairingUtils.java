package com.graphhopper.util;

public class PairingUtils {

    /*
    Szudzik's Function that pairs two numbers to give a unique number.
     */
    public static long pair(int x, int y) {
        long xl = x;
        long yl = y;
        return  x >= y ? xl * xl + xl + yl : yl * yl + xl;
    }

    public static int[] unpair(long z) {
        int b = (int) Math.sqrt(z);
        int a = (int)z - b * b;
        return a < b ? new int[]{a, b} : new int[]{b, a - b};
    }

}
