package com.graphhopper.routing.util.parsers;

import com.graphhopper.reader.ReaderWay;
import com.graphhopper.routing.ev.EncodedValue;
import com.graphhopper.routing.ev.EnumEncodedValue;
import com.graphhopper.routing.ev.TollBicycle;
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
    public void testSimpleTags() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        IntsRef intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollBicycle.MISSING, tollEnc.getEnum(false, intsRef));

        intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:bicycle", "yes");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollBicycle.YES, tollEnc.getEnum(false, intsRef));

        intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:bicycle", "no");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollBicycle.NO, tollEnc.getEnum(false, intsRef));
    }
}