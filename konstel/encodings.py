import base64


def base32(Hash):
	return base64.b32encode(Hash.digest()).decode().lower()


def phoneme_10_5(Hash):
    '''
    Returns phoneme from hashlib Hash instance
    Maps 10 consonants and 5 vowels to hash decimals
    Information density: 2.82 bits per character
    '''
    h = str(int(Hash.hexdigest(), 16))
    vowel_map = dict(zip(map(str, range(10)), 'aeiou'*2))
    consonant_map = dict(zip(map(str, range(10)), 'fhklmprstv'))
    phoneme = ''
    for i, char in enumerate(h):
        if not i % 2:
            phoneme += consonant_map[char]
        else:
            phoneme += vowel_map[char]
    return phoneme


def phoneme_16_4(Hash):
    '''
    Return alternate consonant-vowel phoneme from hash
    Maps 16 consonants and 4 vowels to six bit windows of hash
    Information density: 3 bits per character
    '''
    hash_b2 = bin(int(Hash.hexdigest(), 16))[2:]  # Binary string from hex string
    vowel_map = dict(zip(range(4), 'aiou'))
    consonant_map = dict(zip(range(16), 'bdfghjklmnprstvz'))
    phoneme = ''
    for pos in range(0, len(hash_b2), 6):
        window = hash_b2[pos:pos+6]
        if len(window) < 6:  # Pad with trailing zeros
            window = f'{window:<06}'
        c, v = int(window[:4], 2), int(window[4:], 2)
        phoneme += f'{consonant_map[c]}{vowel_map[v]}'
    return phoneme


def phoneme_bits(seq_hash, result_len=None):
    '''
    Needs integration
    Maps 16 consonants and 4 vowels to hash
    Information density: 3 bits per character
    '''
    import itertools
    from math import log2
    phonemes = [''.join(i) for i in itertools.product('bdfghjklmnprstvz', 'aiou')]
    phoneme_bit_size = log2(len(phonemes))
    assert phoneme_bit_size.is_integer(), 'Must have a power of two number of phonemes'

    hash_b2 = int(seq_hash.hexdigest(), 16)

    if result_len is None:
        result_len = int(8 * seq_hash.digest_size / phoneme_bit_size)
    front_offset = int(8 * seq_hash.digest_size - phoneme_bit_size)
    mask = (2 ** int(phoneme_bit_size) - 1) << front_offset

    word = ''
    for _ in range(result_len):
        phoneme_num = (hash_b2 & mask) >> front_offset
        word += phonemes[phoneme_num]
        hash_b2 = hash_b2 << int(phoneme_bit_size)
        
    return word