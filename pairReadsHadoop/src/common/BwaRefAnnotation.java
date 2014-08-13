// Copyright (C) 2011-2012 CRS4.
//
// This file is part of Seal.
//
// Seal is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// Seal is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.
//
// You should have received a copy of the GNU General Public License along
// with Seal.  If not, see <http://www.gnu.org/licenses/>.


package common;

import java.util.Map;
import java.util.LinkedHashMap;
import java.util.Iterator;
import java.io.LineNumberReader;
import java.io.IOException;
import java.io.Reader;

public class BwaRefAnnotation implements Iterable<BwaRefAnnotation.Contig>
{
	public static class UnknownContigException extends java.lang.RuntimeException {
		private static final long serialVersionUID = 1L;
		public UnknownContigException(String msg) {
			super(msg);
		}
	}

	public static class Contig
	{
		public int id;
		public long start, length;
		public String name;

		public Contig(String name, int id, long start, long length)
		{
			this.id = id;
			this.start = start;
			this.length = length;
			this.name = name;
		}

		public int getId() { return id; }
		public long getStart() { return start; }
		public long getLength() { return length; }
		public String getName() { return name; }
	};

	private enum AnnScannerState {
		NameLine,
		CoordLine
	}

	private Map<String, Contig> contigMap;
	private long referenceLength;

	public BwaRefAnnotation()
	{
		// Use a LinkedHashMap to store the contig information.  This gives us O(1) look-ups
		// by name and preserves the contig order for iteration.  The initial capacity is 30.
		contigMap = new LinkedHashMap<String, Contig>(30);
		referenceLength = -1;
	}

	public BwaRefAnnotation(Reader in) throws IOException
	{
		contigMap = new LinkedHashMap<String, Contig>(30); // initial capacity of 30
		referenceLength = -1;
		this.load(in);
	}

	public void load(Reader in) throws IOException, FormatException
	{
		LineNumberReader input = new LineNumberReader(in);

		String line = null;
		line = input.readLine();
		if (line == null)
			throw new FormatException("Empty annotations file");

		try
		{
			long[] row = scanPosLine(line);

			referenceLength = row[0];
			if (referenceLength <= 0)
				throw new FormatException("Invalid reference length " + referenceLength);
			int nContigs = (int)row[1]; // cast to avoid warning about loss of precision
			if (nContigs <= 0)
				throw new FormatException("Invalid number of contigs " + nContigs);

			AnnScannerState state = AnnScannerState.NameLine;
			int contigCount = 0;
			String lastContigName = null;

			line = input.readLine();
			while (line != null)
			{
				if (line.equals("")) // skip blank lines
					continue;

				if (state == AnnScannerState.NameLine)
				{
					String[] fields = scanNameLine(line);
					contigCount += 1;
					if (contigCount > nContigs)
						throw new FormatException("There are more contigs than expected (first line says we should have " + nContigs + ")");
					lastContigName = fields[1];
					state = AnnScannerState.CoordLine;
				}
				else // state is CoordLine
				{
					long[] fields = scanPosLine(line);
					contigMap.put(lastContigName, new Contig(lastContigName, contigCount, fields[0], fields[1]));
					state = AnnScannerState.NameLine;
				}
				line = input.readLine();
			}
			if (state != AnnScannerState.NameLine)
				throw new FormatException("last entry is incomplete (found the name line but not the coordinates)");
			if (contigCount < nContigs)
				throw new FormatException("Not enough contig records.  Header said we should have " + nContigs + ", but we only found " + contigCount);
		}
		catch (NumberFormatException e) {
			throw new FormatException("Line " + input.getLineNumber() + ": invalid number (" + e.getMessage() + "). Original line: " + line);
		}
		catch (FormatException e) {
			// add line number to message
			throw new FormatException("Line " + input.getLineNumber() + ": " + e.getMessage());
		}
	}

	private long[] scanPosLine(String line) throws NumberFormatException
	{
		String[] fields = line.split("\\s+");
		if (fields.length != 3)
			throw new FormatException("Wrong number of fields (" + fields.length + ").  Expected 3");

		long[] retval = new long[3];
		for (int i = 0; i <= 2; ++i)
		{
			retval[i] = Long.parseLong(fields[i]);
			if (retval[i] < 0)
				throw new NumberFormatException();
		}
		return retval;
	}

	private String[] scanNameLine(String line)
	{
		String[] fields = line.split("\\s+", 3);
		return fields;
	}

	public long getReferenceLength()
	{
		return referenceLength;
	}

	public int getContigId(String name)
	{
		return getContig(name).id;
	}

	public long getAbsCoord(String contig_name, long localCoord)
	{
		Contig contig = getContig(contig_name);
		return contig.start + localCoord;
	}

	public Contig getContig(String name)
	{
		Contig c = contigMap.get(name);
		if (c != null)
			return c;
		else
			throw new UnknownContigException("Unknown contig name '" + name + "'");
	}

	public Iterator<Contig> iterator()
	{
		return contigMap.values().iterator();
	}
}
