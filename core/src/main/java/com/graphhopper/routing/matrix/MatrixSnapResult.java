package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.IntArrayList;
import com.carrotsearch.hppc.IntIntHashMap;
import com.graphhopper.storage.index.Snap;

import java.util.List;


public class MatrixSnapResult {

    List<Snap> snaps;
    private IntArrayList pointsNotFound = new IntArrayList();
    private IntIntHashMap snapIndexes = new IntIntHashMap();

    int size;

    public MatrixSnapResult(List<Snap> snaps, IntIntHashMap snapIndexes,IntArrayList pointsNotFound) {
        this.snaps = snaps;
        this.pointsNotFound = pointsNotFound;
        this.size = snaps.size() + pointsNotFound.size();
        this.snapIndexes = snapIndexes;
    }

    public int size(){
        return size;
    }

    public Snap get(int idx){
        int snapIdx = snapIndexes.get(idx);
        return this.snaps.get(snapIdx);
    }

    public boolean isNotFound(int idx){
        return this.getPointsNotFound().contains(idx);
    }

    public boolean isFound(int idx){
        return !isNotFound(idx);
    }


    public IntArrayList getPointsNotFound() {
        return pointsNotFound;
    }

}
