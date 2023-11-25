import unittest
from unittest.mock import patch, Mock
from vcf_parser import parse_vcf_entry, query_ensembl

class MockVCFEntry:
    def __init__(self, chrom, pos, ref, alt):
        self.CHROM = chrom
        self.POS = pos
        self.REF = ref
        self.ALT = alt

class TestVCFParser(unittest.TestCase):
    def test_parse_vcf_entry(self):
        mock_entry = MockVCFEntry(chrom='1', pos=1000, ref='A', alt=['G'])
        expected_output = {'chrom': '1', 'pos': 1000, 'ref': 'A', 'alt': 'G'}
        self.assertEqual(parse_vcf_entry(mock_entry), expected_output)

    @patch('requests.post')
    def test_query_ensembl(self, mock_post):
        mock_response = Mock()
        mock_response.status_code = 200
        mock_response.json.return_value = {"some_key": "some_value"}
        mock_post.return_value = mock_response

        variant = {'chrom': '1', 'pos': 1000000, 'ref': 'G', 'alt': 'A'}
        response = query_ensembl(variant)
        self.assertEqual(response, {"some_key": "some_value"})

if __name__ == '__main__':
    unittest.main()
