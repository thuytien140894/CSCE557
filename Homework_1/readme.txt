Homework #1: The Vigenere Cipher Decryption
Language/Version: Python 2.7.11

********************************************************************************
COMMAND:
To run the program in the same directory that contains the ciphertext file,
type the command:

python vigenere_decrypter.py ciphertext

********************************************************************************
OUTPUT:
Running the program will print to the console the ciphertext read from the
file, the key used for encryption, and the decrypted text as follows:

Ciphertext:
xkjurowmllpxwznpimbvbqjcnowxpcchhvvfvsllfvxhazityxohulxqojaxelxzxmyjaqfstsrulh
hucdskbxknjqidallpqslluhiaqfpbpcidsvcihwhwewthbtxrljnrsncihuvffuxvoukjljswmaqf
vjwjsdyljogjxdboxajultucpzmpliwmlubzxvoodybafdskxgqfadshxnxehsaruojaqfpfkndhsa
afvulluwtaqfrupwjrszxgpfutjqiynrxnyntwmhcukjfbirzsmehhsjshyonddzzntzmplilrwnmw
mlvuryonthuhabwnvw

Key Used for Encryption:
bdfhj

Decrypted Text:
wheninthecourseofhumaneventsitbecomesnecessaryforonepeopletodissolvethepolitic
albandswhichhaveconnectedthemwithanotherandtoassumeamongthepowersoftheearththe
separateandequalstationtowhichthelawsofnatureandofnaturesgodentitlethemadecent
respecttotheopinionsofmankindrequiresthattheyshoulddeclarethecauseswhichimpelt
hemtotheseparation

********************************************************************************
RESULT:
Inserting appropriate spaces gives the introduction of the Declaration of
Independence:

When in the course of human events it becomes necessary for one people to
dissolve the political bands which have connected them with another and to
assume among the powers of the earth, the separate and equal station to which
the Laws of Nature and of Nature's God entitle them, a decent respect to the
opinions of mankind requires that they should declare the causes which impel
them to the separation.

********************************************************************************
METHOD:
This program attempts to decrypt a message encrypted using the Vigenere Cipher,
which uses a keyword to perform multiple Caesar ciphers based on the letters of
the keyword.  There are three steps to the decryption, finding the key length,
finding the key itself, and decrypting.

I. Find the key length:
First, I shifted the ciphertext one position to its right and line it up with
the original ciphertext.  Then I count the number of letters that coincide
between the two texts.  I repeated this process for 9 times and recorded the
number of coincidences for each n-position shift.

At first, I tried to determine all the possible shifts up to the entire length
of the ciphertext and obtained either 20 or 35 as the possible key lengths.
However, finding the key using the key length 20 resulted in the key with
repeated 5-letter keys.  Then I figured that the key length should be 5 indeed.
Therefore, I decided to limit the shift up to only 9.  The shift of 10 has a
higher number of coincidences and therefore will result in a wrong key if used.
The summary for the number of coincidences for different displacements is as
follows:

displacement:   1   2   3   4   5   6   7   8   9
coincidences:   15  17  13  10  23  5   14  20  5

As observed, the shift of 5 has the highest number of coincidences within the
specified domain and therefore is the best guess for the key length.

II. Find the key:
I used the second method in the textbook to find the key.  This method makes
use of the fact that given the two vectors A(i) and A(j) where A is the vector
representing the frequencies of letters in English, which is:

        A = (0.082, 0.015, 0.028, 0.043, 0.127, 0.022, 0.020, 0.061,
             0.070, 0.002, 0.008, 0.040, 0.024, 0.067, 0.075, 0.019,
             0.001, 0.060, 0.063, 0.091, 0.028, 0.010, 0.023, 0.001,
        		 0.020, 0.001)

And so, A(i) is vector A cyclic right shifted i positions.  For instance, A(3)
is expressed as:

        A(3) = (0.001, 0.020, 0.001, 0.043, 0.127, 0.022, 0.020, 0.061,
                0.070, 0.002, 0.008, 0.040, 0.024, 0.067, 0.075, 0.019,
                0.001, 0.060, 0.063, 0.091, 0.028, 0.010, 0.023)

If we perform a dot product A(i) to itself, the result would be close to 0.066
whereas the dot product of A(i) and A(j) where i != j would range between
0.031 and 0.045.

Based on this fact and the possible key length, I first counted the number of
occurrences of letters at the positions in the ciphertext that match with each
key element.  For example, for the first element of the key, I considered the
letters in the 1st, 6th, 11th, .. positions (in other words, taking the modular
of each position with the key length gives 1).  This vector of occurrences when
divided by the number of letters counted is assumed to be A(i) because it
consists of a random sample of English letters all shifted by the same amount
i.  Our goal is to find A(j) for 0 <= j <= 25 so that the dot product between
A(i) and A(j) is the highest, which implies that i = j and that the original
text at these positions are most likely shifted by j positions. For instance,
he dot products for the first key:

0.0445, 0.0684, 0.0388, 0.0296, 0.0382, 0.0411, 0.0349, 0.0312, 0.0358,
0.0300, 0.0357, 0.0371, 0.0492, 0.0401, 0.0422, 0.0406, 0.0474, 0.0394,
0.0314, 0.0313, 0.0409, 0.0331, 0.0299, 0.0441, 0.0357, 0.0303

The largest value is the second, 0.0684, and so the first key is b.

This process is repeated for as many times as the length of the key.
The key is determined to be (1, 3, 5, 7, 9) = (b, d, f, h, j)

III. Decrypt
Using the found key, I first converted each character in the ciphertext to its
appropriate number, then subtract from it the corresponding key number.  The
decrypted text is described above.
