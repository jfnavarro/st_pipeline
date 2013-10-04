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

import java.util.HashMap;
import java.util.List;
import java.util.ArrayList;
import java.util.Iterator;
import java.nio.ByteBuffer;
import java.io.UnsupportedEncodingException;

/**
 * Abstract readable mapping class with tags.
 */
public abstract class AbstractTaggedMapping
{
	public enum TagDataType {
		Char,
		String,
		Int,
		Float,
		NumArray,
		Bytes;

		public static TagDataType fromSamType(char t)
		{
			switch (t)
			{
				case 'A':
					return TagDataType.Char;
				case 'i':
					return TagDataType.Int;
				case 'f':
					return TagDataType.Float;
				case 'Z':
					return TagDataType.String;
				case 'H':
					return TagDataType.Bytes;
				case 'B':
					return TagDataType.NumArray;
				default:
					throw new IllegalArgumentException("Unknown tag type " + t);
			}
		}
	};

	protected class TagCacheItem {
		private TagDataType type;
		private String value;

		public TagCacheItem(TagDataType type, String value) {
			this.type = type;
			this.value = value;
		}

		public String getValue() { return value; }
		public TagDataType getType() { return type; }
		public String toString() { return "(" + type + "," + value + ")"; }
	}

	////////////////////////////////////////////////
	// variables
	////////////////////////////////////////////////
	protected HashMap<String, TagCacheItem> tagCache;

	public AbstractTaggedMapping()
	{
		tagCache = new HashMap<String, TagCacheItem>();
	}

	////////////////////////////////////////////////
	// methods
	////////////////////////////////////////////////
	abstract public String getName();
	abstract public int getFlag();
	abstract public String getContig() throws IllegalStateException;
	abstract public int get5Position() throws IllegalStateException;
	abstract public int getMapQ();

	abstract public List<AlignOp> getAlignment() throws IllegalStateException;

	/**
	 * This mapping's DNA sequence.
	 * The ASCII representation of the base sequence is contained in the ByteBuffer,
	 * starting at buffer.position() and ending at buffer.limit() (exclusive).
	 * The buffer is mark()ed at the start of the sequence.
	 */
	abstract public ByteBuffer getSequence();

	/**
	 * Return this mapping's sequence as a string.
	 *
	 * A new string is created from the internal ByteBuffer and returned.
	 */
	public String getSequenceString()
	{
		return byteBufferToString(getSequence());
	}

	/**
	 * This mapping's DNA sequence's base quality scores.
	 * The ASCII Phred+33 (Sanger encoding) representation of the base quality
	 * scores sequence is contained in the ByteBuffer,
	 * starting at buffer.position() and ending at buffer.limit() (exclusive).
	 * The buffer is mark()ed at the start of the sequence.
	 */
	abstract public ByteBuffer getBaseQualities();

	/**
	 * Return this mapping's base qualities as a string.
	 *
	 * A new string is created from the internal ByteBuffer and returned.
	 */
	public String getBaseQualitiesString()
	{
		return byteBufferToString(getBaseQualities());
	}

	/**
	 * Get the sequence's length.
	 */
	abstract public int getLength();

	/**
	 * Returns a string representation of this read mapping.
	 * The string looks something like a partial SAM format.  It contains the following
	 * tab-separated fields:
	 * name, flag, contig, 5' position, sequence.
	 */
	public String toString()
	{
		StringBuilder builder = new StringBuilder(1000);
		builder.
			append(getName()).append("\t").
			append(getFlag()).append("\t");

		if (isMapped())
		{
			builder.
			  append(getContig()).append("\t").
			  append(get5Position()).append("\t");
		}
		else
			builder.append("*\t*\t");

		ByteBuffer seq = getSequence();

		try {
			builder.append( new String(seq.array(), seq.position(), seq.limit() - seq.position(), "US-ASCII") );
		}
		catch (UnsupportedEncodingException e) {
			throw new RuntimeException("??? US-ASCII charset not supported! " + e);
		}

		return builder.toString();
	}

	//////////////////////// tag methods ////////////////////////


	/**
	 * Retrieve a TagCacheItem for the named tag.
	 * Checks the cache and retrieves the item from the cache if it's there.  If
	 * it's not, it calls getTagText to retrieve it, scans it creating a
	 * TagCacheItem which it caches and returns to the caller.
	 * @return a TagCacheItem for the given tag name.  The returned object is the one in this object's cache.
	 * @exception NoSuchFieldException Tag not found in mapping record.
	 * @exception FormatException Invalid SAM tag syntax.
	 */
	abstract protected TagCacheItem makeTagItem(String name) throws NoSuchFieldException;

	/**
	 * Get a tag item, going through the cache.
	 *
	 * This method checks the cache to see if the relevant tag item is already present.
	 * If it isn't, it calls makeTagItem(name) to fetch it, if possibile, and
	 * then saves it new item in the cache before returning it to the caller.
	 */
	protected TagCacheItem getTagItem(String name) throws NoSuchFieldException
	{
		TagCacheItem item = tagCache.get(name);
		if (item == null)
		{
			item = makeTagItem(name);
			tagCache.put(name, item);
		}
		return item;
	}

	public String getTag(String name) throws NoSuchFieldException
	{
		return getTagItem(name).getValue();
	}

	public boolean hasTag(String name)
	{
		try {
			getTagItem(name);
			return true;
		}
		catch (NoSuchFieldException e) {
			return false;
		}
	}

	public int getIntTag(String name) throws NoSuchFieldException
	{
		TagCacheItem item = getTagItem(name);
		if (item.getType() == TagDataType.Int)
			return Integer.parseInt(item.getValue());
		else
			throw new NumberFormatException("item " + item + " is not of integer type");
	}

	public double getDoubleTag(String name) throws NoSuchFieldException
	{
		TagCacheItem item = getTagItem(name);
		if (item.getType() == TagDataType.Float || item.getType() == TagDataType.Int)
			return Double.parseDouble(item.getValue());
		else
			throw new NumberFormatException("item " + item + " is not of double type");
	}

	//////////////////////// flag methods ////////////////////////

	public boolean isPaired() {
		return AlignFlags.Paired.is(getFlag());
	}

	public boolean isProperlyPaired() {
		return AlignFlags.ProperlyPaired.is(getFlag());
	}

	public boolean isMapped() {
		return AlignFlags.Unmapped.isNot(getFlag());
	}

	public boolean isUnmapped() {
		return AlignFlags.Unmapped.is(getFlag());
	}

	public boolean isMateMapped() {
		return AlignFlags.MateUnmapped.isNot(getFlag());
	}

	public boolean isMateUnmapped() {
		return AlignFlags.MateUnmapped.is(getFlag());
	}

	public boolean isOnReverse() {
		return AlignFlags.OnReverse.is(getFlag());
	}

	public boolean isMateOnReverse() {
		return AlignFlags.MateOnReverse.is(getFlag());
	}

	public boolean isRead1() {
		return AlignFlags.Read1.is(getFlag());
	}

	public boolean isRead2() {
		return AlignFlags.Read2.is(getFlag());
	}

	public boolean isSecondaryAlign() {
		return AlignFlags.SecondaryAlignment.is(getFlag());
	}

	public boolean isFailedQC() {
		return AlignFlags.FailedQC.is(getFlag());
	}

	public boolean isDuplicate() {
		return AlignFlags.Duplicate.is(getFlag());
	}

	//////////////////////// template/insert size methods ////////////////////////

	/**
	 * Check whether the template length is available.
	 * If it returns true the method getTemplateLength() will not throw.
	 */
	abstract public boolean isTemplateLengthAvailable();

	/**
	 * Get the observed template length.
	 * The value returned from this method is closely related to the "TLEN: signed
	 * observed template length" field in the SAM spec, but <b>slighly different</b>:
	 * the value is always &gt; 0.  When the template length is unavailable this method
	 * throws an IllegalStateException.  This can happen when:
	 * <ul>
	 *   <li>the mapping is for an unpaired read;</li>
	 *   <li>either this read or the mate are unmapped;</li>
	 *   <li>read and mate are mapped to different contigs.</li>
	 * </ul>
	 *
	 * @exception IllegalStateException Template length is unavailable.
	 */
	abstract public int getTemplateLength() throws IllegalStateException;

	//////////////////////// utility methods ////////////////////////
	protected static String byteBufferToString(ByteBuffer buffer)
	{
		try {
			return new String(buffer.array(), buffer.position(), buffer.limit() - buffer.position(), "US-ASCII");
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("Can't get US-ASCII encoding!");
		}
	}

	// These below should probably be in another class.

	/**
	 * Calculate the reference coordinate for each matched position in this mapping.
	 * This read must not be unmapped.
	 *
	 * Writes getLength() integers into dest.  Non-matched positions (AlignOp != Match)
	 * are set to -1, while the rest are set to the reference position to which
	 * the base at the corresponding position was matched.
	 *
	 */
	public void calculateReferenceCoordinates(ArrayList<Integer> dest) throws IllegalStateException
	{
		if (isUnmapped())
			throw new IllegalStateException("Can't calculate reference coordinates for an unmapped read");

		dest.clear();
		List<AlignOp> alignment = getAlignment();

		int refPos = get5Position();

		for (AlignOp op: alignment)
		{
			int opLength = op.getLen();
			if (op.getType() == AlignOp.Type.Match)
			{
				int endPos = refPos + opLength;
				for ( ; refPos < endPos; ++refPos)
					dest.add(refPos);
			}
			else if (op.getType() == AlignOp.Type.Insert || op.getType() == AlignOp.Type.SoftClip)
			{
				for (int i = op.getLen(); i > 0; --i)
					dest.add(-1);
			}
			else if (op.getType() == AlignOp.Type.Delete)
			{
				refPos += op.getLen();
			}
			else
				throw new RuntimeException("Found " + op.getType() + " alignment operation.  Sorry, but I don't know how to deal with this.");
		}

		if (dest.size() != getLength())
			throw new RuntimeException("Inconsistency?  Alignment in " + AlignOp.cigarStr(alignment) + " doesn't exactly cover the entire read (got " + dest.size() + " positions for " + getLength() + " bases)");
	}

	/**
	 * Calculate which read positions match or don't match their corresponding reference position.
	 * Calculation is based on the read's alignment and MD tag.
	 * This object must be mapped.  The dest array will be cleared and then filled with one value
	 * per base in the read.  True indicates that the base matches the reference; false indicates
	 * a mismatch; null indicates that the base doesn't have a matching reference base (it's an
	 * insertion or it has been clipped).
	 *
	 * @exception IllegalStateException This mapping is unmapped.
	 * @exception RuntimeException Unexpected errors such as bad CIGAR or MD tags.
	 */
	public void calculateReferenceMatches(ArrayList<Boolean> dest) throws IllegalStateException
	{
		if (isUnmapped())
			throw new IllegalStateException("Can't calculate reference coordinates for an unmapped read");

		dest.clear();
		dest.ensureCapacity(getLength());

		List<AlignOp> alignment = getAlignment();
		if (alignment.isEmpty())
			throw new RuntimeException("No alignment for read " + this);

		try
		{
			List<MdOp> md = MdOp.scanMdTag(getTag("MD"));
			if (md.isEmpty())
				throw new RuntimeException("no MD operators extracted from tag! (tag: " + getTag("MD") + ")");

			Iterator<MdOp> mdIt = md.iterator();

			MdOp mdOp = mdIt.next();
			int mdOpConsumed = 0; // the number of positions within the current mdOp that have been consumed.

			// we iterate through the alignment.  We really only care about operations which
			// "consume" bases from our read:  insertions, soft clippings and matches.
			//
			// In the case of inserts and soft clips we don't have a corresponding reference base,
			// so we leave the match value as null.
			//
			// In the case of a cigar Match, we consult the MD string.  We advance along the MD matches
			// and mismatches inserting corresponding true or false values.
			for (AlignOp alignOp: alignment)
			{
				if (alignOp.getType() == AlignOp.Type.Match)
				{
					int positionsToCover = alignOp.getLen();
					while (positionsToCover > 0 && mdOp != null)
					{
						if (mdOp.getType() == MdOp.Type.Delete)
							mdOp = mdIt.next(); // skip it
						else
						{
							// must be a match or a mismatch
							boolean match = mdOp.getType() == MdOp.Type.Match;
							int consumed = Math.min(mdOp.getLen() - mdOpConsumed, positionsToCover);
							for (int i = consumed; i > 0; --i)
								dest.add(match);
							positionsToCover -= consumed;
							mdOpConsumed += consumed;
							if (mdOpConsumed >= mdOp.getLen()) // operator consumed.  Advance to next
							{
								mdOpConsumed = 0;
								if (mdIt.hasNext())
									mdOp = mdIt.next();
								else
									mdOp = null;
							}
						}
					}
					if (positionsToCover > 0 && mdOp == null)
						throw new RuntimeException("BUG or bad data?? Found more read positions than was covered by the MD tag. CIGAR: " + AlignOp.cigarStr(alignment) + "; MD: " + getTag("MD") + "; read: " + this.toString());
				}
				else if (alignOp.getType() == AlignOp.Type.Insert || alignOp.getType() == AlignOp.Type.SoftClip)
				{
					for (int i = alignOp.getLen(); i > 0; --i)
						dest.add(null);
				}
				// else the op doesn't affect the read
			}
		} catch (NoSuchFieldException e) {
			throw new RuntimeException(e.getMessage());
		}
	}

}
