// Copyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Martin C. Frith

#include "Alignment.hh"
#include "GeneticCode.hh"
#include "LastEvaluer.hh"
#include "MultiSequence.hh"
#include "Alphabet.hh"
#include <iomanip>
#include <algorithm>
#include <cassert>
#include <cstdio>  // sprintf

// make C++ tolerable:
#define CI(type) std::vector<type>::const_iterator

using namespace cbrc;

// write x - y as a signed integer
static void writeSignedDifference( size_t x, size_t y, std::ostream& os ){
  if( x >= y )  os << x - y;
  else          os << '-' << y - x;
}

void Alignment::write( const MultiSequence& seq1, const MultiSequence& seq2,
		       char strand, bool isTranslated, const Alphabet& alph,
		       const LastEvaluer& evaluer, int format,
		       std::ostream& os, const AlignmentExtras& extras ) const{
  assert( !blocks.empty() );
  if( format == 0 ) writeTab( seq1, seq2, strand, isTranslated,
			      evaluer, os, extras );
  else              writeMaf( seq1, seq2, strand, isTranslated, alph,
			      evaluer, os, extras );
}

static char* writeTaggedItems( const LastEvaluer& evaluer, double queryLength,
			       int score, double fullScore,
			       char separator, char* out ){
  if( evaluer.isGood() ){
    double epa = evaluer.evaluePerArea( score );
    double area = evaluer.area( score, queryLength );
    *out++ = separator;
    out += std::sprintf( out, "EG2=%.2g", 1e18 * epa );
    *out++ = separator;
    out += std::sprintf( out, "E=%.2g", area * epa );
  }
  if( fullScore > 0 ){
    *out++ = separator;
    out += std::sprintf( out, "fullScore=%.3g", fullScore );
  }
  *out++ = '\n';
  return out;
}

void Alignment::writeTab( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, bool isTranslated,
			  const LastEvaluer& evaluer, std::ostream& os,
			  const AlignmentExtras& extras ) const{
  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t w1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(w1);

  size_t size2 = seq2.finishedSize();
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  size_t seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  size_t s1 = seq1.seqLen(w1);
  size_t s2 = seq2.seqLen(w2);

  os << score << '\t';

  os << seq1.seqName(w1) << '\t'
     << alnBeg1 - seqStart1 << '\t'
     << alnEnd1 - alnBeg1 << '\t'
     << '+' << '\t'
     << s1 << '\t';

  os << seq2.seqName(w2) << '\t'
     << alnBeg2 - seqStart2 << '\t'
     << alnEnd2 - alnBeg2 << '\t'
     << strand << '\t'
     << s2 << '\t';

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;
      if( j->size ) os << ',';
      size_t gapBeg1 = j->end1();
      size_t gapEnd1 = i->beg1();
      writeSignedDifference( gapEnd1, gapBeg1, os );  // allow -1 frameshift
      os << ':';
      size_t gapBeg2 = aaToDna( j->end2(), frameSize2 );
      size_t gapEnd2 = aaToDna( i->beg2(), frameSize2 );
      writeSignedDifference( gapEnd2, gapBeg2, os );  // allow -1 frameshift
      if( i->size ) os << ',';
    }
    if( i->size ) os << i->size;
  }

  double fullScore = extras.fullScore;
  char line[256];
  char* end = writeTaggedItems( evaluer, s2, score, fullScore, '\t', line );
  os.write( line, end - line );
}

static int numDigits( size_t x ){
  int n = 0;
  do{
    ++n;
    x /= 10;
  }while(x);
  return n;
}

// Printing with either C++ streams or sprintf can be noticeably slow.
// So the next 3 functions are used instead.

static char* sprintLeft( char* dest, const char* src, int width ){
  char* end = dest + width;
  while( *src ) *dest++ = *src++;
  while( dest < end ) *dest++ = ' ';
  *dest++ = ' ';
  return dest;
}

static char* sprintSize( char* dest, size_t size, int width ){
  char* end = dest + width;
  char* beg = end;

  do{
    --beg;
    *beg = '0' + size % 10;
    size /= 10;
  }while( size );

  while( dest < beg ) *dest++ = ' ';

  *end++ = ' ';
  return end;
}

static char* sprintChar( char* dest, char c ){
  *dest++ = c;
  *dest++ = ' ';
  return dest;
}

void Alignment::writeMaf( const MultiSequence& seq1, const MultiSequence& seq2,
			  char strand, bool isTranslated, const Alphabet& alph,
			  const LastEvaluer& evaluer, std::ostream& os,
			  const AlignmentExtras& extras ) const{
  double fullScore = extras.fullScore;
  const std::vector<uchar>& columnAmbiguityCodes = extras.columnAmbiguityCodes;
  const std::vector<double>& expectedCounts = extras.expectedCounts;

  size_t alnBeg1 = beg1();
  size_t alnEnd1 = end1();
  size_t w1 = seq1.whichSequence(alnBeg1);
  size_t seqStart1 = seq1.seqBeg(w1);

  size_t size2 = seq2.finishedSize();
  size_t frameSize2 = isTranslated ? (size2 / 3) : 0;
  size_t alnBeg2 = aaToDna( beg2(), frameSize2 );
  size_t alnEnd2 = aaToDna( end2(), frameSize2 );
  size_t w2 = seq2.whichSequence( strand == '+' ? alnBeg2 : size2 - alnBeg2 );
  size_t seqStart2 = strand == '+' ? seq2.seqBeg(w2) : size2 - seq2.seqEnd(w2);

  const std::string n1 = seq1.seqName(w1);
  const std::string n2 = seq2.seqName(w2);
  size_t b1 = alnBeg1 - seqStart1;
  size_t b2 = alnBeg2 - seqStart2;
  size_t r1 = alnEnd1 - alnBeg1;
  size_t r2 = alnEnd2 - alnBeg2;
  size_t s1 = seq1.seqLen(w1);
  size_t s2 = seq2.seqLen(w2);

  const int nw = std::max( n1.size(), n2.size() );
  const int bw = std::max( numDigits(b1), numDigits(b2) );
  const int rw = std::max( numDigits(r1), numDigits(r2) );
  const int sw = std::max( numDigits(s1), numDigits(s2) );

  size_t headLen = 2 + nw + 1 + bw + 1 + rw + 3 + sw + 1;
  size_t lineLen = headLen + numColumns( frameSize2 ) + 1;
  std::vector<char> lineVector( lineLen );
  char* line = &lineVector[0];
  char* tail = line + headLen;
  line[ lineLen - 1 ] = '\n';
  char* dest;

  char aLine[256];
  dest = sprintChar( aLine, 'a' );
  dest += std::sprintf( dest, "score=%d", score );
  dest = writeTaggedItems( evaluer, s2, score, fullScore, ' ', dest );
  os.write( aLine, dest - aLine );

  dest = sprintChar( line, 's' );
  dest = sprintLeft( dest, n1.c_str(), nw );
  dest = sprintSize( dest, b1, bw );
  dest = sprintSize( dest, r1, rw );
  dest = sprintChar( dest, '+' );
  dest = sprintSize( dest, s1, sw );
  writeTopSeq( seq1.seqReader(), alph, 0, frameSize2, dest );
  os.write( line, lineLen );

  size_t qualsPerBase1 = seq1.qualsPerLetter();
  if( qualsPerBase1 ){
    dest = sprintChar( line, 'q' );
    dest += nw + 1;
    std::fill( dest, tail, ' ' );
    writeTopSeq( seq1.qualityReader(), alph, qualsPerBase1, frameSize2, tail );
    os.write( line, lineLen );
  }

  dest = sprintChar( line, 's' );
  dest = sprintLeft( dest, n2.c_str(), nw );
  dest = sprintSize( dest, b2, bw );
  dest = sprintSize( dest, r2, rw );
  dest = sprintChar( dest, strand );
  dest = sprintSize( dest, s2, sw );
  writeBotSeq( seq2.seqReader(), alph, 0, frameSize2, dest );
  os.write( line, lineLen );

  size_t qualsPerBase2 = seq2.qualsPerLetter();
  if( qualsPerBase2 && !isTranslated ){
    // for translated alignment: don't write untranslated quality data
    dest = sprintChar( line, 'q' );
    dest += nw + 1;
    std::fill( dest, tail, ' ' );
    writeBotSeq( seq2.qualityReader(), alph, qualsPerBase2, frameSize2, tail );
    os.write( line, lineLen );
  }

  if( columnAmbiguityCodes.size() > 0 ){
    dest = sprintChar( line, 'p' );
    std::fill( dest, tail, ' ' );
    copy( columnAmbiguityCodes.begin(), columnAmbiguityCodes.end(), tail );
    os.write( line, lineLen );
  }

  if( expectedCounts.size() > 0 ){
    os << 'c';
    for( unsigned i = 0; i < expectedCounts.size(); ++i )
      os << ' ' << expectedCounts[i];
    os << '\n';
  }

  os << '\n';  // blank line afterwards
}

size_t Alignment::numColumns( size_t frameSize ) const{
  size_t num = 0;

  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // length of unaligned chunk of top sequence (gaps in bottom sequence):
      num += i->beg1() - j->end1();

      // length of unaligned chunk of bottom sequence (gaps in top sequence):
      size_t gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) ++num;
      num += gap2;
    }

    num += i->size;  // length of aligned chunk
  }

  return num;
}

static char* writeGaps( char* dest, size_t num ){
  char* end = dest + num;
  while( dest < end ){
    *dest++ = '-';
  }
  return dest;
}

static char* writeQuals( const uchar* qualities, size_t beg, size_t end,
			 size_t qualsPerBase, char* dest ){
  for( size_t i = beg; i < end; ++i ){
    const uchar* q = qualities + i * qualsPerBase;
    *dest++ = *std::max_element( q, q + qualsPerBase );
  }
  return dest;
}

static char* writeSeq( const uchar* seq, size_t beg, size_t end,
		       const Alphabet& alph, size_t qualsPerBase, char* dest ){
  return qualsPerBase ? writeQuals( seq, beg, end, qualsPerBase, dest )
    :                   alph.rtCopy( seq + beg, seq + end, dest );
}

char* Alignment::writeTopSeq( const uchar* seq, const Alphabet& alph,
			      size_t qualsPerBase, size_t frameSize,
			      char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // append unaligned chunk of top sequence:
      dest = writeSeq( seq, j->end1(), i->beg1(), alph, qualsPerBase, dest );

      // append gaps for unaligned chunk of bottom sequence:
      size_t gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 ) *dest++ = '-';
      dest = writeGaps( dest, gap2 );
    }

    // append aligned chunk of top sequence:
    dest = writeSeq( seq, i->beg1(), i->end1(), alph, qualsPerBase, dest);
  }

  return dest;
}

char* Alignment::writeBotSeq( const uchar* seq, const Alphabet& alph,
			      size_t qualsPerBase, size_t frameSize,
			      char* dest ) const{
  for( CI(SegmentPair) i = blocks.begin(); i < blocks.end(); ++i ){
    if( i > blocks.begin() ){  // between each pair of aligned blocks:
      CI(SegmentPair) j = i - 1;

      // append gaps for unaligned chunk of top sequence:
      dest = writeGaps( dest, i->beg1() - j->end1() );

      //append unaligned chunk of bottom sequence:
      size_t gap2, frameshift2;
      sizeAndFrameshift( j->end2(), i->beg2(), frameSize, gap2, frameshift2 );
      if( frameshift2 == 1 ) *dest++ = '\\';
      if( frameshift2 == 2 ) *dest++ = '/';
      size_t chunkBeg2 = i->beg2() - gap2;
      dest = writeSeq( seq, chunkBeg2, i->beg2(), alph, qualsPerBase, dest );
    }

    // append aligned chunk of bottom sequence:
    dest = writeSeq( seq, i->beg2(), i->end2(), alph, qualsPerBase, dest );
  }

  return dest;
}
