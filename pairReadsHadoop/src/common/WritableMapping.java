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

import java.util.List;
import java.util.Set;
import java.nio.ByteBuffer;
import java.io.UnsupportedEncodingException;

/**
 * Simple in-memory mapping
 */
public class WritableMapping extends AbstractTaggedMapping
{
	////////////////////////////////////////////////
	// variables
	////////////////////////////////////////////////
	protected String        name;
	protected int           flag       = 0x4; // unmapped bit set.  All the rest at 0
	protected String        contig;
	protected int           pos        = 0;
	protected int           mapq       = 255;
	protected List<AlignOp> alignment;
	protected ByteBuffer    sequence;
	protected ByteBuffer    quality;
	protected int           insertSize = 0;

	public WritableMapping() { }

	/**
	 * Construct a simple mapping from ASCII-encoded sequence and base quality scores.
	 */
	public WritableMapping(String name, String seq, String qual)
	{
		if (seq.length() != qual.length())
			throw new IllegalArgumentException(String.format("sequence and quality strings have different lengths! (%d and %d)", seq.length(), qual.length()));

		this.name = name;
		setSequence(seq);
		setBaseQualities(qual);
	}

	public WritableMapping(String name, ByteBuffer seq, ByteBuffer qual)
	{
		this.name = name;
		this.sequence = seq;
		this.quality = qual;
	}
	
	////////////////////////////////////////////////
	// methods
	////////////////////////////////////////////////
	public String getName() { return name; }
	public void setName(String name) { this.name = name; }

	public int getFlag() { return flag; }

	public void setFlag(int flag)
	{
		this.flag = flag;

		if (isUnmapped())
			mapq = 0;

		if (isUnmapped() || !isPaired() || isMateUnmapped())
			insertSize = 0;
	}

	public String getContig() throws IllegalStateException
	{
		if (isUnmapped())
			throw new IllegalStateException();
		return contig;
	}

	public void setContig(String contig) { this.contig = contig; }

	public int get5Position() throws IllegalStateException
	{
		if (isUnmapped())
			throw new IllegalStateException();
	 	return pos;
	}

	public void set5Position(int p) { pos = p; }

	public int getMapQ()
	{
		if (isUnmapped())
			throw new IllegalStateException();
	 	return mapq;
	}

	public void setMapQ(int q)
	{
		if (q < 0 || q > 255)
			throw new IllegalArgumentException("mapq must be between 0 and 255 (got " + q + ")");
		mapq = q;
	}

	public List<AlignOp> getAlignment() throws IllegalStateException
	{
		if (isUnmapped())
			throw new IllegalStateException();

		return alignment;
	}
	
	public void setAlignment(List<AlignOp> align) { alignment = align; }

	public ByteBuffer getSequence() { return sequence; }

	/**
	 * Set buf as this mapping's sequence.
	 *
	 * A reference to buf is kept by this object.
	 */
	public void setSequence(ByteBuffer buf) { sequence = buf; }

	/**
	 * Set this mapping's sequence to the contents of seq.
	 *
	 * The byte array seq is copied to this object's internal ByteBuffer.
	 */
	public void setSequence(byte[] seq, int offset, int length)
	{
		sequence = ByteBuffer.allocate(length);
		sequence.put(seq, offset, length).rewind().mark();
	}

	public void setSequence(byte[] seq)
	{
		setSequence(seq, 0, seq.length);
	}

	/**
	 * Set this mapping's sequence to the contents of seq.
	 *
	 * A new ByteBuffer is created from the contents of the seq string and
	 * it is set as this mapping's sequence.
	 */
	public void setSequence(String seq)
	{
		try {
			sequence = ByteBuffer.wrap(seq.getBytes("US-ASCII"));
			sequence.mark();
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("Can't get US-ASCII charset!");
		}
	}

	public ByteBuffer getBaseQualities() { return quality; }
	
	/**
	 * Set buf as this mapping's base qualities.
	 *
	 * A reference to buf is kept by this object.
	 */
	public void setBaseQualities(ByteBuffer buf) { quality = buf; }

	/**
	 * Set this mapping's base qualities to the contents of seq.
	 *
	 * The byte array seq is copied to this object's internal ByteBuffer.
	 */
	public void setBaseQualities(byte[] seq, int offset, int length)
	{
		quality = ByteBuffer.allocate(length);
		quality.put(seq, offset, length).rewind().mark();
	}

	/**
	 * Set this mapping's base qualities to the contents of seq.
	 *
	 * The byte array seq is copied to this object's internal ByteBuffer.
	 */
	public void setBaseQualities(byte[] seq)
	{
		setBaseQualities(seq, 0, seq.length);
	}

	/**
	 * Set this mapping's base qualities to the contents of seq.
	 *
	 * A new ByteBuffer is created from the contents of the seq string and
	 * it is set as this mapping's base qualities.
	 */
	public void setBaseQualities(String seq)
	{
		try {
			quality = ByteBuffer.wrap(seq.getBytes("US-ASCII"));
			quality.mark();
		} catch (UnsupportedEncodingException e) {
			throw new RuntimeException("Can't get US-ASCII charset!");
		}
	}

	public int getLength()
	{
		return sequence.limit() - sequence.position();
	}

	/**
	 * Clear this mapping.
	 *
	 * Returns it to the state of a new mapping allocated with the default constructor.
	 */
	public void clear()
	{
		name       = null;
		flag       = 0x4;
		contig     = null;
		pos        = 0;
		mapq       = 255;
		alignment  = null;
		sequence   = null;
		quality    = null;
		insertSize = 0;
		tagCache.clear();
	}

	protected TagCacheItem makeTagItem(String name) throws NoSuchFieldException
	{
		throw new NoSuchFieldException("no field named " + name);
	}

	public void setTag(String name, TagDataType type, String value)
	{
		tagCache.put(name, new TagCacheItem(type, value));
	}

	public Set<String> getTagNames()
	{
		return tagCache.keySet();
	}

	//////////////////////// template/insert size methods ////////////////////////

	public void setTemplateLength(int isize)
	{
		if (isize != 0)
		{
			if (isUnmapped() || !isPaired() || isMateUnmapped())
				throw new IllegalStateException("Setting non-zero template size for read in incompatible state");
		}
		insertSize = isize;
	}

	public boolean isTemplateLengthAvailable() { return insertSize != 0; }

	public int getTemplateLength() throws IllegalStateException
	{
		if (!isTemplateLengthAvailable())
			throw new IllegalStateException();
		return insertSize;
	}

	//////////////////////// flag set methods ////////////////////////

	public void setIsPaired(boolean v) {
		setFlag(v, AlignFlags.Paired);
	}

	public void setIsProperlyPaired(boolean v) {
		setFlag(v, AlignFlags.ProperlyPaired);
	}

	public void setIsMapped(boolean v) {
		setFlag(!v, AlignFlags.Unmapped);
	}

	public void setIsUnmapped(boolean v) {
		setFlag(v, AlignFlags.Unmapped);
	}

	public void setIsMateMapped(boolean v) {
		setFlag(!v, AlignFlags.MateUnmapped);
	}

	public void setIsMateUnmapped(boolean v) {
		setFlag(v, AlignFlags.MateUnmapped);
	}

	public void setIsOnReverse(boolean v) {
		setFlag(v, AlignFlags.OnReverse);
	}

	public void setIsMateOnReverse(boolean v) {
		setFlag(v, AlignFlags.MateOnReverse);
	}

	public void setIsRead1(boolean v) {
		setFlag(v, AlignFlags.Read1);
	}

	public void setIsRead2(boolean v) {
		setFlag(v, AlignFlags.Read2);
	}

	public void setIsSecondaryAlign(boolean v) {
		setFlag(v, AlignFlags.SecondaryAlignment);
	}

	public void setIsFailedQC(boolean v) {
		setFlag(v, AlignFlags.FailedQC);
	}

	public void setIsDuplicate(boolean v) {
		setFlag(v, AlignFlags.Duplicate);
	}

	protected void setFlag(boolean newValue, AlignFlags f)
	{
		if (newValue)
			flag = f.set(flag);
		else
			flag = f.clear(flag);
	}
}
