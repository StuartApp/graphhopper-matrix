package com.graphhopper.routing.util.parsers;

import com.graphhopper.reader.ReaderWay;
import com.graphhopper.routing.ev.EncodedValue;
import com.graphhopper.routing.ev.EnumEncodedValue;
import com.graphhopper.routing.ev.TollMotorcycle;
import com.graphhopper.storage.IntsRef;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class OSMTollMotorcycleParserTest {
    private EnumEncodedValue<TollMotorcycle> tollEnc;
    private OSMTollMotorcycleParser parser;

    @BeforeEach
    public void setUp() {
        tollEnc = new EnumEncodedValue<>(TollMotorcycle.KEY, TollMotorcycle.class);
        tollEnc.init(new EncodedValue.InitializerConfig());
        parser = new OSMTollMotorcycleParser(tollEnc);
    }

    @Test
    public void testSimpleTags() {
        ReaderWay readerWay = new ReaderWay(1);
        IntsRef relFlags = new IntsRef(2);
        IntsRef intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollMotorcycle.MISSING, tollEnc.getEnum(false, intsRef));

        intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:motorcycle", "yes");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollMotorcycle.YES, tollEnc.getEnum(false, intsRef));

        intsRef = new IntsRef(1);
        readerWay.setTag("highway", "primary");
        readerWay.setTag("toll:motorcycle", "no");
        parser.handleWayTags(intsRef, readerWay, relFlags);
        assertEquals(TollMotorcycle.NO, tollEnc.getEnum(false, intsRef));
    }
}