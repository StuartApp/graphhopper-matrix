package com.graphhopper.speeds;

import java.util.Objects;

public class SpeedKmByHour {

    private double speed;

    public SpeedKmByHour(double speed) {
        this.speed = speed;
    }

    public double getSpeed() {
        return speed;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        SpeedKmByHour that = (SpeedKmByHour) o;
        return Double.compare(speed, that.speed) == 0;
    }

    @Override
    public int hashCode() {
        return Objects.hash(speed);
    }

    @Override
    public String toString() {
        return "SpeedKmByHour{" +
                "speed=" + speed +
                '}';
    }
}
