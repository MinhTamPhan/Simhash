package com.vng.zing.simhash;
import com.google.common.collect.Lists;
import it.unimi.dsi.fastutil.longs.LongOpenHashSet;
import it.unimi.dsi.fastutil.longs.LongSet;

import java.nio.CharBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import javafx.util.Pair;

/**
 * 27-09-2018
 * MinhTam(Tampm2)
 */
public class SimHash {
    public static final int  HASH_SIZE          = 64;
    public static final long HASH_RANGE         = 2 ^ HASH_SIZE;
    public static MurmurHash hasher             = new MurmurHash();
    public static long weakBit = 0;
    /**
     * use short cuts to obtains a speed optimized simhash calculation
     *
     * @param s
     *          input string
     * @return 64 bit simhash of input string
     */

    private static final int FIXED_CGRAM_LENGTH = 4;

    public static long computeOptimizedSimHashForString(String s) {
        return computeOptimizedSimHashForString(CharBuffer.wrap(s));
    }

    public static long computeOptimizedSimHashForString(CharBuffer s) {

        LongSet shingles = new LongOpenHashSet(Math.min(s.length(), 1000000));

        int length = s.length();

        long timeStart = System.currentTimeMillis();
        for (int i = 0; i < length - FIXED_CGRAM_LENGTH + 1; i++) {
            // extract an ngram

            long shingle = s.charAt(i);
            shingle <<= 16;
            shingle |= s.charAt(i + 1);
            shingle <<= 16;
            shingle |= s.charAt(i + 2);
            shingle <<= 16;
            shingle |= s.charAt(i + 3);

            shingles.add(shingle);
        }
        long timeEnd = System.currentTimeMillis();
        // System.out.println("NGram Production Took:" + (timeEnd-timeStart));

        int v[] = new int[HASH_SIZE];
        byte longAsBytes[] = new byte[8];

        for (long shingle : shingles) {

            longAsBytes[0] = (byte) (shingle >> 56);
            longAsBytes[1] = (byte) (shingle >> 48);
            longAsBytes[2] = (byte) (shingle >> 40);
            longAsBytes[3] = (byte) (shingle >> 32);
            longAsBytes[4] = (byte) (shingle >> 24);
            longAsBytes[5] = (byte) (shingle >> 16);
            longAsBytes[6] = (byte) (shingle >> 8);
            longAsBytes[7] = (byte) (shingle);

            long longHash = FPGenerator.std64.fp(longAsBytes, 0, 8);
            for (int i = 0; i < HASH_SIZE; ++i) {
                boolean bitSet = ((longHash >> i) & 1L) == 1L;
                v[i] += (bitSet) ? 1 : -1;
            }
        }

        long simhash = 0;
        for (int i = 0; i < HASH_SIZE; ++i) {
            if (v[i] > 0) {
                simhash |= (1L << i);
            }
        }
        
        weakBit = simhash >> 56;
        return simhash;
    }

    public static int hammingDistance(long hash1, long hash2) {
        long bits = hash1 ^ hash2;
        int count = 0;
        while (bits != 0) {
            bits &= bits - 1;
            ++count;
        }
        return count;
    }
    
    public static float ProbabilisticNearDuplicate(long hash1, long hash2) {
        return 1.0f - (hammingDistance(hash1, hash2) * 1.0f / 64);
    }

    public static long rotate(long hashValue) {
        return (hashValue << 1) | (hashValue >>> -1);
    }
    
    public static Pair<Long, Integer[]> computeSimHashForString(CharBuffer s) {

        LongSet shingles = new LongOpenHashSet(Math.min(s.length(), 1000000));

        int length = s.length();
        for (int i = 0; i < length - FIXED_CGRAM_LENGTH + 1; i++) {
            // extract an ngram

            long shingle = s.charAt(i);
            shingle <<= 16;
            shingle |= s.charAt(i + 1);
            shingle <<= 16;
            shingle |= s.charAt(i + 2);
            shingle <<= 16;
            shingle |= s.charAt(i + 3);

            shingles.add(shingle);
        }
        
        int v[] = new int[HASH_SIZE];
        byte longAsBytes[] = new byte[8];

        for (long shingle : shingles) {

            longAsBytes[0] = (byte) (shingle >> 56);
            longAsBytes[1] = (byte) (shingle >> 48);
            longAsBytes[2] = (byte) (shingle >> 40);
            longAsBytes[3] = (byte) (shingle >> 32);
            longAsBytes[4] = (byte) (shingle >> 24);
            longAsBytes[5] = (byte) (shingle >> 16);
            longAsBytes[6] = (byte) (shingle >> 8);
            longAsBytes[7] = (byte) (shingle);

            long longHash = FPGenerator.std64.fp(longAsBytes, 0, 8);
            for (int i = 0; i < HASH_SIZE; ++i) {
                boolean bitSet = ((longHash >> i) & 1L) == 1L;
                v[i] += (bitSet) ? 1 : -1;
            }
        }

        long simhash = 0;
        Integer[] bit = new Integer[HASH_SIZE];
        for (int i = 0; i < HASH_SIZE; ++i) {
            if (v[i] > 0) {
                simhash |= (1L << i);
            }
            bit[i] = v[i];
        }
       
        
        return new Pair<Long, Integer[]>(simhash, bit);
    }

    public static List<Long> getHearderBit(List<Integer> frequencyBit, long simHash){
        List<List<Integer>> weakBits = buildWeakbits(frequencyBit);
        List<Long> result = new ArrayList<>();
        result.add(simHash >> 56);
        for (List<Integer> weakBit : weakBits) {
            long tmp = simHash ^ (1 << (64 - weakBit.get(0)));
            tmp = tmp ^ (1 << (64 - weakBit.get(1)));
            tmp = tmp ^ (1 << (64 - weakBit.get(2)));
            if(!result.contains(tmp)) {
                result.add(tmp);
            }
        }
        return result;
    }
    
    private static  List<List<Integer>> buildWeakbits(List<Integer> frequencyBit){
        List<Long> result = new ArrayList<>();
        List<Integer> sortList = Lists.newArrayList(frequencyBit);
        sortList.sort((o1, o2) -> {
            return Integer.compare(Math.abs(o1), Math.abs(o2));
        });
        List<Integer> weaks = new ArrayList<>();
        weaks.add(frequencyBit.indexOf(sortList.get(0)));
        int j = 0;
        for (int i = 1; i < 7; i++) {
            j = frequencyBit.indexOf(sortList.get(i));
            //case if frequencyBit[i] == frequencyBit[i + 1]
            while (weaks.contains(j)) {
                ++j;
            }
            weaks.add(j);
        }
        List<List<Integer>> asList = 
                Arrays.asList(Arrays.asList(weaks.get(0), weaks.get(1), weaks.get(2)),
                            Arrays.asList(weaks.get(0), weaks.get(1), weaks.get(3)),
                            Arrays.asList(weaks.get(0), weaks.get(1), weaks.get(4)),
                            Arrays.asList(weaks.get(0), weaks.get(2), weaks.get(3)),
                            Arrays.asList(weaks.get(0), weaks.get(1), weaks.get(5)),
                            Arrays.asList(weaks.get(0), weaks.get(2), weaks.get(4)),
                            Arrays.asList(weaks.get(1), weaks.get(2), weaks.get(3)),
                            Arrays.asList(weaks.get(0), weaks.get(1), weaks.get(6)),
                            Arrays.asList(weaks.get(0), weaks.get(2), weaks.get(5)),
                            Arrays.asList(weaks.get(0), weaks.get(3), weaks.get(4)),
                            Arrays.asList(weaks.get(1), weaks.get(2), weaks.get(4)));
        
//        for (List<Integer> list : asList) {
//            long tmp = 0;
//            tmp |= (1 << list.get(0));
//            tmp |= (1 << list.get(1));
//            tmp |= (1 << list.get(2));
//            result.add(tmp);
//        }
        return asList;
    }

}
