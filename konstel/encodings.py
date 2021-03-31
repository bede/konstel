import base64
import itertools
from math import ceil, log


def base32(string_hash):
    '''Returns lowercased and padding-stripped base32 encoding of a hashlib Hash'''
    return base64.b32encode(string_hash.digest()).decode().lower().rstrip('=')


def cbase32(string_hash):
    '''Returns lowercased Crockford's base32 encoding of a hashlib Hash'''
    hash_b32 = base32(string_hash)
    base32_symbols = 'abcdefghijklmnopqrstuvwxyz234567'
    cbase32_symbols = '0123456789abcdefghjkmnpqrstvwxyz'
    base32_to_cbase32 = dict(zip(base32_symbols, cbase32_symbols))
    return ''.join([base32_to_cbase32[c] for c in hash_b32])


def phonemes_10_5(string_hash):
    '''
    Returns word comprising consonant-vowel phonemes from hashlib Hash
    Maps 10 consonants and 5 vowels to base10 encoded hash
    2.8 bits per character
    '''
    hash_b10 = str(int(string_hash.hexdigest(), 16))
    vowel_map = dict(zip(map(str, range(10)), 'aeiou'*2))
    consonant_map = dict(zip(map(str, range(10)), 'fhklmprstv'))
    word = ''
    for i, char in enumerate(hash_b10):
        if not i % 2:
            word += consonant_map[char]
        else:
            word += vowel_map[char]
    return word


def phonemes_16_4(string_hash):
    '''
    Returns word comprising consonant-vowel phonemes from hashlib Hash
    Maps 16 consonants and 4 vowels to six bit windows of hash
    3 bits per character
    '''
    hash_b2 = bin(int(string_hash.hexdigest(), 16))[2:]  # Binary string from hex string
    vowel_map = dict(zip(range(4), 'aiou'))
    consonant_map = dict(zip(range(16), 'bdfghjklmnprstvz'))
    word = ''
    if len(hash_b2) % 6 != 0:
        # Pad with leading zeros
        hash_b2 = hash_b2.zfill(6 * ceil(len(hash_b2) / 6))
    for pos in range(0, len(hash_b2), 6):
        window = hash_b2[pos:pos+6]
        c, v = int(window[:4], 2), int(window[4:], 2)
        word += f'{consonant_map[c]}{vowel_map[v]}'
    return word


def phonemes_16_4_bits_old(string_hash, result_len=None):
    '''
    Returns word comprising consonant-vowel phonemes from hashlib Hash
    Maps 16 consonants and 4 vowels to six bit windows of hash
    3 bits per character

    Uses bit-wise logic to calculate which phonemes are which.
    '''
    from math import log2
    phonemes = [''.join(i) for i in itertools.product('bdfghjklmnprstvz', 'aiou')]
    phoneme_bit_size = log2(len(phonemes))
    assert phoneme_bit_size.is_integer(), 'Must have a power of two number of phonemes'

    hash_b2 = int(string_hash.hexdigest(), 16)

    if result_len is None:
        result_len = int(8 * string_hash.digest_size / phoneme_bit_size)
    front_offset = int(phoneme_bit_size * result_len - phoneme_bit_size)
    mask = (2 ** int(phoneme_bit_size) - 1) << front_offset

    word = ''
    for _ in range(result_len):
        phoneme_num = (hash_b2 & mask) >> front_offset
        word += phonemes[phoneme_num]
        hash_b2 = hash_b2 << int(phoneme_bit_size)
    return word


def phonemes_16_4_bits(string_hash, result_len=None, phonemes=None):
    '''
    Returns word comprising `string_hash` represented by `phonemes` of max
    number `result_len`.

    This function transforms the `string_hash` into a number of base
    `len(phonemes)` and then returns that number printed out with each digit
    represented as a phoneme. This has the advantage, relative to other phoneme
    schemes, of allowing phoneme lists of non-base-2 length.
    '''
    if phonemes is None:
        phonemes = [''.join(i) for i in itertools.product('bdfghjklmnprstvz', 'aiou')]
    hash_b2 = int(string_hash.hexdigest(), 16)
    if result_len is None:
        result_len = ceil(log(hash_b2, len(phonemes)))

    word = []
    for _ in range(ceil(log(hash_b2, len(phonemes)))):
        phoneme_num, hash_b2 = hash_b2 % len(phonemes), hash_b2 // len(phonemes)
        word.append(phonemes[phoneme_num])
    return ''.join(word[:-result_len-1:-1])
