"""Independent, standard-library-only checker for integral homology class receipts."""
import argparse
from hashlib import sha256
import json
from math import gcd
from pathlib import Path

_KIND = "integral-homology-class/1"
_FIELDS = {
    "noncycle": {"row", "value"},
    "boundary": {"lift"},
    "infinite": {"functional"},
    "torsion": {"order", "lift", "witnesses"},
}


def _vector(value, width):
    return type(value) in (list, tuple) and len(value) == width and all(
        type(entry) is int for entry in value
    )


def _pair(vector, other):
    return sum(a * b for a, b in zip(vector, other, strict=True))


def _normalize(ranks, boundaries):
    if type(ranks) not in (tuple, list) or any(
        type(rank) is not int or rank < 0 for rank in ranks
    ):
        raise ValueError("invalid chain ranks")
    if type(boundaries) is not dict or any(
        type(n) is not int or not 1 <= n < len(ranks) for n in boundaries
    ):
        raise ValueError("invalid boundary degrees")

    normalized = {}
    for n in range(1, len(ranks)):
        height, width = ranks[n - 1], ranks[n]
        if n not in boundaries:
            normalized[n] = [[0] * width for _ in range(height)]
            continue
        matrix = boundaries[n]
        if type(matrix) not in (list, tuple) or len(matrix) != height or any(
            not _vector(row, width) for row in matrix
        ):
            raise ValueError(f"invalid matrix in degree {n}")
        normalized[n] = [list(row) for row in matrix]

    for n in range(2, len(ranks)):
        if any(
            _pair(row, column)
            for column in zip(*normalized[n], strict=True)
            for row in normalized[n - 1]
        ):
            raise ValueError(f"nonzero differential composition in degree {n}")
    return normalized


def verify_class_receipt(ranks, boundaries, degree, chain, receipt):
    """Return True for a valid receipt. Raise ValueError for invalid base input.

    Ranks, row-major differentials, degree and chain are supplied independently
    of the receipt. Missing differentials represent zero maps.
    """
    maps = _normalize(ranks, boundaries)
    if type(degree) is not int or not 0 <= degree < len(ranks):
        raise ValueError("invalid chain degree")
    width = ranks[degree]
    if not _vector(chain, width):
        raise ValueError("invalid chain coordinates")

    if type(receipt) is not dict:
        return False
    kind = receipt.get("kind")
    if type(kind) is not str or kind != _KIND:
        return False
    status = receipt.get("status")
    if type(status) is not str or status not in _FIELDS:
        return False
    if receipt.keys() != {"kind", "input_sha256", "status"} | _FIELDS[status]:
        return False
    encoded = json.dumps([ranks, maps, degree, chain], sort_keys=True, separators=(",", ":"))
    digest = receipt["input_sha256"]
    if type(digest) is not str or digest != sha256(encoded.encode("ascii")).hexdigest():
        return False

    outgoing = [_pair(row, chain) for row in maps.get(degree, ())]
    if status == "noncycle":
        row, value = receipt["row"], receipt["value"]
        return (
            type(row) is int and 0 <= row < len(outgoing)
            and type(value) is int and value != 0 and outgoing[row] == value
        )
    if any(outgoing):
        return False

    incoming = maps.get(degree + 1, [()] * width)
    columns = tuple(zip(*incoming, strict=True))
    if status == "infinite":
        ell = receipt["functional"]
        return (
            _vector(ell, width) and _pair(ell, chain) != 0
            and all(_pair(ell, column) == 0 for column in columns)
        )

    order = 1 if status == "boundary" else receipt["order"]
    if type(order) is not int or (status == "torsion" and order < 2):
        return False
    lift = receipt["lift"]
    lift_width = ranks[degree + 1] if degree + 1 < len(ranks) else 0
    if not _vector(lift, lift_width) or any(
        _pair(row, lift) != order * coordinate
        for row, coordinate in zip(incoming, chain, strict=True)
    ):
        return False
    if status == "boundary":
        return True

    witnesses = receipt["witnesses"]
    if type(witnesses) is not list or not witnesses:
        return False
    # Remove from the claimed order the factors already forced by earlier rows.
    remaining = order
    covered = 1
    for witness in witnesses:
        if type(witness) is not dict or witness.keys() != {"modulus", "functional"}:
            return False
        modulus, ell = witness["modulus"], witness["functional"]
        if type(modulus) is not int or modulus < 2 or not _vector(ell, width):
            return False
        if any(_pair(ell, column) % modulus for column in columns):
            return False
        image_order = modulus // gcd(modulus, _pair(ell, chain))
        new_factor = image_order // gcd(image_order, covered)
        if remaining % new_factor:
            return False
        remaining //= new_factor
        covered *= new_factor
    return remaining == 1


def _unique_object(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate JSON key {key!r}")
        result[key] = value
    return result


def _reject_number(token):
    raise ValueError(f"expected integer JSON numbers, got {token!r}")


def _read_json(path):
    return json.loads(
        path.read_text(encoding="utf-8"), object_pairs_hook=_unique_object,
        parse_float=_reject_number, parse_constant=_reject_number,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="trusted ranks, boundaries, degree and chain JSON")
    parser.add_argument("receipt", type=Path, help="untrusted integral class receipt JSON")
    args = parser.parse_args()

    try:
        problem = _read_json(args.input)
        if type(problem) is not dict or problem.keys() != {"ranks", "boundaries", "degree", "chain"}:
            raise ValueError("invalid input fields")
        raw_maps = problem["boundaries"]
        if type(raw_maps) is not dict:
            raise ValueError("boundaries must be a JSON object")
        maps = {}
        for key, value in raw_maps.items():
            degree = int(key)
            if str(degree) != key:
                raise ValueError("boundary degrees must use canonical decimal keys")
            maps[degree] = value
        valid = verify_class_receipt(
            problem["ranks"], maps, problem["degree"], problem["chain"], _read_json(args.receipt),
        )
    except (OSError, ValueError, RecursionError) as error:
        parser.exit(2, f"Invalid input: {error}\n")
    print("valid receipt" if valid else "invalid receipt")
    return 0 if valid else 1


if __name__ == "__main__":
    raise SystemExit(main())
