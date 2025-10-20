from rest_framework import permissions

class IsAdminOrReadOnly(permissions.BasePermission):
    """
    Public bisa baca (GET), admin bisa ubah data.
    """

    def has_permission(self, request, view):
        # kalau request GET, HEAD, atau OPTIONS, boleh untuk semua
        if request.method in permissions.SAFE_METHODS:
            return True

        # kalau bukan GET (berarti POST/PUT/DELETE), hanya admin boleh
        user = request.user
        return user.is_authenticated and (user.is_staff or getattr(user, 'role', '') == 'admin')