package com.graphhopper;

import com.graphhopper.config.Profile;
import com.graphhopper.routing.ev.*;
import com.graphhopper.routing.util.AllEdgesIterator;
import com.graphhopper.speeds.SpeedKmByHour;
import com.graphhopper.speeds.WaySpeedsProvider;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.util.List;
import java.util.Optional;
import java.util.concurrent.TimeUnit;
import java.util.stream.Collectors;

public class GraphHopperCustomSpeeds extends GraphHopper {

    private static final Logger logger = LoggerFactory.getLogger(GraphHopperCustomSpeeds.class);

    private WaySpeedsProvider speedsProvider;

    public GraphHopperCustomSpeeds(WaySpeedsProvider speedsProvider) {

        this.speedsProvider = speedsProvider;
    }

    @Override
    protected void importPublicTransit() {

        logger.info("Start Custom Speeds Provider");
        long startTime = System.nanoTime();

        List<String> vehicles = this.getProfiles().stream().map(Profile::getVehicle).distinct().collect(Collectors.toList());

        List<DecimalEncodedValue> speedEncoders = vehicles.stream().map(
                vehicle -> this.getEncodingManager().getDecimalEncodedValue(VehicleSpeed.key(vehicle))).collect(Collectors.toList());
        setSpeedsFor(speedEncoders);

        long endTime = System.nanoTime();
        long elapsed = (endTime - startTime);
        long durationInMillis = TimeUnit.NANOSECONDS.toMillis(elapsed);
        logger.info("End Custom Speeds Provider : Graph processed in " + durationInMillis);
    }

    private Optional<Double> speedFor(int osmWayId, RoadClass roadClass) {

        Optional<Double> waySpeed = speedsProvider.speedForWay(osmWayId).map(SpeedKmByHour::getSpeed);

        if (waySpeed.isPresent()) {
            return waySpeed;
        } else {
            return speedsProvider.speedForRoadClass(roadClass).map(SpeedKmByHour::getSpeed);
        }
    }

    private void setSpeedsFor(List<DecimalEncodedValue> speedEncoders) {
        IntEncodedValue OSMWayIDEncoder = this.getEncodingManager().getIntEncodedValue(OSMWayID.KEY);
        DecimalEncodedValue maxSpeedEncoder = this.getEncodingManager().getDecimalEncodedValue(MaxSpeed.KEY);
        EnumEncodedValue<RoadClass> roadClassEncoder = this.getEncodingManager().getEnumEncodedValue(RoadClass.KEY, RoadClass.class);

        AllEdgesIterator iter = this.getBaseGraph().getAllEdges();
        while (iter.next()) {

            int edge = iter.getEdge();

            //Normal Direction
            double maxSpeed = iter.get(maxSpeedEncoder);
            int osmWayId = iter.get(OSMWayIDEncoder);
            RoadClass roadClass = iter.get(roadClassEncoder);
            Optional<Double> maybeCustomSpeed = speedFor(osmWayId, roadClass);

            if (maybeCustomSpeed.isPresent()) {
                double customSpeed = maybeCustomSpeed.get();
                double speedToSet;
                if (customSpeed > maxSpeed) {
                    logger.warn("OsmId " + osmWayId + ": Custom Speed (" + customSpeed + ") > MaxSpeed ( " + maxSpeed + ") for edge reverse " + edge + ", so max_speed will be used");
                    speedToSet = maxSpeed;
                } else {
                    speedToSet = customSpeed;
                }
                speedEncoders.forEach(encoder -> {
                    double maxEncoderValue = encoder.getMaxStorableDecimal();
                    double minEncoderValue = encoder.getMinStorableDecimal();
                    if(speedToSet <= maxEncoderValue && speedToSet >= minEncoderValue){
                        double speed = iter.get(encoder);
                        logger.debug("Replace " + speed + " with " + speedToSet + " for edge " + edge);
                        iter.set(encoder, speedToSet);
                    }else{
                        logger.warn("OsmId " + osmWayId + ": Invalid Custom Speed (" + speedToSet + ")  for encoder " + encoder.getName() + " ( " + minEncoderValue + " - " + maxEncoderValue+ ") for edge " + edge + ",so it will be ignored");
                    }
                });

            }

            //Reverse Direction
            int osmWayIdReverse = iter.getReverse(OSMWayIDEncoder);
            double maxSpeedReverse = iter.getReverse(maxSpeedEncoder);
            RoadClass roadClassReverse = iter.getReverse(roadClassEncoder);
            Optional<Double> maybeCustomSpeedReverse = speedFor(osmWayIdReverse, roadClassReverse);

            if (maybeCustomSpeedReverse.isPresent()) {
                double customSpeedReverse = maybeCustomSpeedReverse.get();
                double speedToSet;
                if (customSpeedReverse > maxSpeedReverse) {
                    logger.warn("OsmId " + osmWayId + ": Custom Speed (" + customSpeedReverse + ") > MaxSpeed ( " + maxSpeedReverse + ") for edge reverse " + edge + ", so max_speed will be used");
                    speedToSet = maxSpeedReverse;
                } else {
                    speedToSet = customSpeedReverse;
                }
                speedEncoders.forEach(encoder -> {
                    if (encoder.isStoreTwoDirections()) {
                        double maxEncoderValue = encoder.getMaxStorableDecimal();
                        double minEncoderValue = encoder.getMinStorableDecimal();
                        if(speedToSet <= maxEncoderValue && speedToSet >= minEncoderValue){
                            double speed = iter.getReverse(encoder);
                            logger.debug("Replace " + speed + " with " + speedToSet + " for edge reverse " + edge);
                            iter.setReverse(encoder, speedToSet);
                        }else{
                            //logger.warn("OsmId " + osmWayId + ": Invalid Custom Speed (" + speedToSet + ") for encoder " + encoder.getName() + " ( " + minEncoderValue + " - " + maxEncoderValue+ ") for edge reverse " + edge + ",so it will be ignored");
                        }
                    }
                });
            }

        }
    }
}
