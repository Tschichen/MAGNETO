
	def bk_multiple_alignment(self, matching_b, save_all, random_pivot=False, score=False, user_list_score=False, nodes_edges_both=False):
		mod_prod = self.modular_product(matching_b)
		checked = []
		clique_list = []
		clique = nx.Graph()
		if random_pivot:
			self.BK_random(mod_prod, clique, list(mod_prod.nodes), checked, clique_list)
		else:
			self.BK(mod_prod, clique, list(mod_prod.nodes), checked, clique_list)
		if not score:
			max_clique = self.find_largest_clique(clique_list)
		else:
			max_clique = keep.score_from_BK(clique_list, user_list_score, nodes_edges_both):
		self = self.add_nodes_to_alignment(matching_b, max_clique)
		if save_all:
			try:
				io.GraphIO.write_graphML_file(self, save_all)
			except:
				print("Could not write to " + save_all + " !")

		return self
		
	def bk_multiple_alignment_clique(self, matching_b, save_all, all_clique, score=False, user_list_score=None, nodes_edges_both=None):
		mod_prod = self.modular_product(matching_b)
		# p_without_forbidden = self.modular_product_valid_matchings(mod_prod, score_list)
		checked = []
		clique_list = []
		clique = nx.Graph()
		matching_name_a = list(self.network.nodes(data=True))[0][1]['Graph']
		matching_name_b = list(matching_b.network.nodes(data=True))[0][1]['Graph']

		clique1 = all_clique.clique_dict[matching_name_a]
		clique2 = all_clique.clique_dict[matching_name_b]
		clique = all_clique.get_nodes(clique1, clique2, mod_prod)
		clique = all_clique.get_edges(clique1, clique2, clique)
		candidates = all_clique.filter_candidates_mult(clique1, clique2, clique, mod_prod)

		clique_list.append(clique)
		self.BK(mod_prod, clique, candidates, checked, clique_list)
		if not score:
			max_clique = self.find_largest_clique(clique_list)
		else:
			max_clique = keep.score_from_BK(clique_list, user_list_score, nodes_edges_both)
		self = self.add_nodes_to_alignment(matching_b, max_clique)
		if save_all:
			try:
				io.GraphIO.write_graphML_file(self, save_all)
			except:
				print("Could not write to " + save_all + " !")

		return self
		
			# alignment_order = [('c', 'b', 'cb'), ('cb', 'a', 'cba'), ('d', 'cba', 'dcba')]
	# Update für Matching Objekte
	def MGA_bk(self, matching_object_dict, list_alignment_order, save_all=None, pre_clique=False, random_pivot=False, score=False, user_list_score=None, nodes_edges_both=None):
		for i in range(len(list_alignment_order)):
			if save_all:
				save_all_as = save_all + str(i + 1) + ".graphml"
			else:
				save_all_as = None
			matching_name_a = list_alignment_order[i][0]
			matching_name_b = list_alignment_order[i][1]
			print("alignment step: " + str(matching_name_a) + ", " + str(matching_name_b))
			# a und b aus dem dict suchen
			matching_a = matching_object_dict[matching_name_a]
			matching_b = matching_object_dict[matching_name_b]
			if not random_pivot:
				if not pre_clique:
					if not score:
						matching_ab = matching_a.bk_multiple_alignment(matching_b, save_all_as)
						matching_object_dict[list_alignment_order[i][2]] = matching_ab
					else:
						matching_ab = matching_a.bk_multiple_alignment(matching_b, save_all_as, score, user_list_score, nodes_edges_both)
						matching_object_dict[list_alignment_order[i][2]] = matching_ab
				else:
					if not score:
						matching_ab = matching_a.bk_multiple_alignment_clique(matching_b, save_all_as, pre_clique)
						matching_object_dict[list_alignment_order[i][2]] = matching_ab
					else:
						matching_ab = matching_a.bk_multiple_alignment_clique(matching_b, save_all_as, pre_clique, score, user_list_score, nodes_edges_both)
						matching_object_dict[list_alignment_order[i][2]] = matching_ab
			else:
				if not score:
					matching_ab = matching_a.bk_multiple_alignment(matching_b, save_all_as, random_pivot)
					matching_object_dict[list_alignment_order[i][2]] = matching_ab
				else:
					matching_ab = matching_a.bk_multiple_alignment(matching_b, save_all_as, random_pivot, score, user_list_score, nodes_edges_both)
					matching_object_dict[list_alignment_order[i][2]] = matching_ab