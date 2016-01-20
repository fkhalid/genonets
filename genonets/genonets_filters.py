
"""
    genonets_filters
    ~~~~~~~~~~~~~~~~

    Contains filters used throughout the package.

    :author: Fahad Khalid
    :license: MIT, see LICENSE for more details.
"""

class WriterFilter :
	@staticmethod
	def gmlAttribsToIgnore(level) :
		if level == "network" :
			attrs = [
				"Evolvability_targets",
				"SqrEpi_list",
				"diameterPath_list",
				"Squares_list"
			]
		elif level == "vertex" :
			attrs = [
				"label",					
				"pathsToSummit",
				"VtxToSqrs"
			]

		return attrs

	@staticmethod
	def netAttribsToIgnore() :
		return	[
					"name",
					"SqrEpi_list",
					"diameterPath_list",
					"Squares_list",
					"Summit"
				]

	@staticmethod
	def seqAttribsToIgnore() :
		return	[
					"sequences",
					"label",
					"escores",
					"pathsToSummit",
					"VtxToSqrs"
				]
