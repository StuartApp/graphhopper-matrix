package com.graphhopper.routing.matrix;

import com.carrotsearch.hppc.IntArrayList;
import com.graphhopper.storage.index.Snap;

import java.util.List;


public class MatrixSnapResult {

    List<Snap> snaps;
    private IntArrayList pointsNotFound = new IntArrayList();

    int size;

    public MatrixSnapResult(List<Snap> snaps, IntArrayList pointsNotFound) {
        this.snaps = snaps;
        this.pointsNotFound = pointsNotFound;
        this.size = snaps.size() + pointsNotFound.size();
    }

    public int size(){
        return size;
    }

    public Snap get(int idx){
        return this.snaps.get(idx);
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
