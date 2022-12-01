/*
 *  Licensed to GraphHopper GmbH under one or more contributor
 *  license agreements. See the NOTICE file distributed with this work for
 *  additional information regarding copyright ownership.
 *
 *  GraphHopper GmbH licenses this file to you under the Apache License,
 *  Version 2.0 (the "License"); you may not use this file except in
 *  compliance with the License. You may obtain a copy of the License at
 *
 *       http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */
package com.graphhopper.routing.matrix;

import java.util.Objects;

public class Bucket {

    public double weight;
    public double distance;
    public long time;
    public int idx;



    public Bucket(double weight, long time, double distance, int idx) {
        this.weight = weight;
        this.time = time;
        this.distance = distance;
        this.idx = idx;
    }

    @Override
    public String toString() {
        return "BucketEntry{" +
                "weight=" + weight +
                ", distance=" + distance +
                ", time=" + time +
                ", idx=" + idx +
                '}';
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Bucket that = (Bucket) o;
        return idx == that.idx;
    }

    @Override
    public int hashCode() {
        return Objects.hash(idx);
    }
}