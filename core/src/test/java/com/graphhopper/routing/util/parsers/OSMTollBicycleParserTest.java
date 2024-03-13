package com.graphhopper.routing.util.parsers;

import com.graphhopper.reader.ReaderWay;
import com.graphhopper.routing.ev.*;
import com.graphhopper.storage.IntsRef;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class OSMTollBicycleParserTest {
    private EnumEncodedValue<TollBicycle> tollEnc;
    private OSMTollBicycleParser parser;

    @BeforeEach
    public void setUp() {
        tollEnc = new EnumEncodedValue<>(TollBicycle.KEY, TollBicycle.class);
        tollEnc.init(new EncodedValue.InitializerConfig());
        parser = new OSMTollBicycleParser(tollEnc);
    }

    @Test
    public void testTollBicycleYes() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:bicycle", "yes");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollBicycle.YES, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }

    @Test
    public void testTollBicycleNo() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;

        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:bicycle", "no");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollBicycle.NO, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }

    @Test
    public void testTollBicycleMissing() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        EdgeIntAccess edgeIntAccess = new ArrayEdgeIntAccess(1);
        int edgeId = 0;

        readerWay.setTag("highway", "primary");
        parser.handleWayTags(edgeId, edgeIntAccess, readerWay, relFlags);
        assertEquals(TollBicycle.MISSING, tollEnc.getEnum(false, edgeId, edgeIntAccess));
    }
}
