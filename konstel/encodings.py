import base64
import binascii


def base32(Hash):
	return base64.b32encode(Hash.digest()).decode().lower()


def phoneme_b10(Hash):
    h_b10 = str(int(binascii.hexlify(Hash.digest()), 16))
    vowel_map = dict(zip(map(str, range(10)), 'aeiou'*2))
    consonant_map = dict(zip(map(str, range(10)), 'fhklmprstv'))
    phoneme = ''
    for char_i, char in enumerate(h_b10):
        if not char_i % 2:
            phoneme += consonant_map[char]
        else:
            phoneme += vowel_map[char]
    return phoneme